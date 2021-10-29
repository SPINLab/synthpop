library('fitdistrplus')
library('FNN')
library("ggplot2")
library('mixtools')

set.seed(47)

# Author: Xander Wilcke
# Email: w.x.wilcke@vu.nl
#
# R script to estimate the parameters of the distributions underlying the data.
# To ensure that no information from the datasets is leaked, or can be restored,
# we use Monte Carlo sampling, the adding of noise, and the default assumption
# that the distribution can be represented as a Gaussian mixture model.

computeGMMs <- function(df_data) {
	# remember the original data types
	df_types <- data.frame(matrix(ncol=6, nrow=0,
								  dimnames=list(NULL, c("name", "type", "num_values",
														"num_distinct", "terms_min", "terms_max"))))
	for (i in 1:length(colnames(df_data))) {
		terms_min <- NA
		terms_max <- NA
		if (class(df_data[,i])=="factor") {
			levels_trimmed <- sapply(levels(df_data[,i]), trimws)  # trim whitespace
			num_terms <- lengths(regmatches(levels_trimmed, gregexpr("\\s", levels_trimmed)))

			terms_min <- min(num_terms) + 1
			terms_max <- max(num_terms) + 1
		}

		df_types[i,] <- c(colnames(df_data)[i], class(df_data[,i]), sum(!is.na(df_data[,i])),
						  length(unique(df_data[,i])), terms_min, terms_max)
	}
	class(df_types$num_values) <- "numeric"
	class(df_types$num_distinct) <- "numeric"
	
	# output table
	num_added <- 0
	df_out <- data.frame(matrix(ncol=11, nrow=0,
		    				    dimnames=list(NULL, c("name", "type", "num_values", "num_distinct",
													  "terms_min", "terms_max", "dist", "min", "max",
													  "mu", "sigma"))))

	list_out <- list()  # to store plots

	# convert to numerical values
	df_data[] <- lapply(df_data, as.numeric)

	for (c in colnames(df_data)) {
		message(paste("- estimating distribution on column", c))
		if (identical(c, "RINPERSOON") || identical(c, "RINPERSOONS")) {
			message("-- skipped identifiers")
			next
		}


		c_info <- df_types[df_types$name==c,]
		if (c_info$num_values < 256 || c_info$num_distinct <= 1) {
			# skip column if too few values
			# we cope with few distinct values by adding noise
			message("-- skipped because too few (distinct) values")
			next
		}

		tryCatch({
			num_samples <- 10e4
			samples <- sample(df_data[[c]][!is.na(df_data[[c]])], size=num_samples, replace=TRUE)

			# add a bit of Gaussian noise to obscure actual data
			samples <- addNoise(samples)

			# estimate if data is distributed uniformly, else assume normal
			samples_norm <- as.vector(scale(samples))
			fit_unif <- fitdist(samples_norm, "unif")$aic
			fit_norm <- fitdist(samples_norm, "norm")$aic

			if (fit_unif < fit_norm) {
				message("-- guessed uniform; computing statistics...")
				out <- approxUniform(samples)
			} else {
				message("-- guessed Gaussian; computing statistics...")
				out <- approxNormal(samples, samples_norm)
			}

			df_values <- out[[1]]
			plt <- out[[2]] + ggtitle(c)
			list_out[[c]] <- plt

			from <- num_added
			to <- from + nrow(df_values) - 1
			j <- 1
			for (i in from:to) {
				df_out[i+1,] <- c(c_info, df_values[j,])

				j <- j + 1
				num_added <- num_added + 1
			}

		}, error = function(e) { message(paste("--- failed with error:", (e$message))) }
		)
	}

	return (list(df_out, list_out))
}

addNoise <- function(samples, multiplier=0.05) {
	samples <- samples + rnorm(length(samples), mean(samples) * multiplier, sd(samples) * multiplier)

	return (samples)
}

approxUniform <- function(samples) {
	df_out <- data.frame(matrix(ncol=5, nrow=0,
							  dimnames=list(NULL, c("dist", "min", "max", "mu", "sigma"))))

	df_out[1,] <- c("unif", min(samples), max(samples), NA, NA)

	plt <- create_hist(samples)

	return (list(df_out, plt))
}

isMixture <- function(samples, tolerance=0.01) {
	# estimate if data is a mixture model by computing the KL divergence
	# between the actual distribution and that of a distribution that
	# is expected if it is not a mixture model.
	tryCatch({
		n <- length(samples)
		mu <- mean(samples)
		sigma <- sd(samples)

		div_expected <- mean(KL.divergence(rnorm(n, mean=mu, sd=sigma),
										   rnorm(n, mean=mu, sd=sigma)))
		div_actual <- mean(KL.divergence(samples,
										 rnorm(n, mean=mu, sd=sigma)))

		if (isTRUE(all.equal.numeric(div_expected, div_actual, tolerance=tolerance))) {
			return (FALSE)
		}

		return (TRUE)
	}, error = function(e) { message(paste("----- failed with error:", (e$message))) }
	)

	return (FALSE)  # safe fallback
}

estimateNumClusters <- function(samples, num_clusters.max=8, sq_penalty_weight=1e-4) {
	n <- length(samples)
	m <- floor(n/2)

	message("---- estimating clusters")
	# df to store results
	KLdiv <- data.frame(k=c(2:num_clusters.max))
	KLdiv[, "divergence"] <- NA

	# Monte Carlo approximation to KL divergence for multivariate
	# Gaussians to estimate optimal value of k:
	for (k in 2:num_clusters.max) {
		message(paste("----- estimating number of clusters: trying", k))
		tryCatch({
			samples <- sample(samples)  # shuffle

			# compute two GMMs on different parts of the data
			gmm_a <- normalmixEM(samples[1:m], k=k, maxrestarts=5, fast=TRUE)
			gmm_b <- normalmixEM(samples[m:n], k=k, maxrestarts=5, fast=TRUE)

			# draw new samples from found distributions
			samples_a <- rnorm(n, mean=gmm_a$mu, sd=gmm_a$sigma)
			samples_b <- rnorm(n, mean=gmm_b$mu, sd=gmm_b$sigma)

			# compute KL divergence between distributions 
			KLdiv[k-1, 2] <- abs(mean(KL.divergence(samples_a, samples_b)))
		}, error = function(e) { message(paste("--- failed with error:", (e$message))) }
		)
	}

	# add penalty to higher k
	KLdiv$divergence <- KLdiv$divergence + (KLdiv$k ** 2) * sq_penalty_weight

	# which k yielded the lowest value of KL divergence
	k_opt <- KLdiv[which.min(KLdiv$divergence),]$k

	print(KLdiv)
	
	message(paste("---- optimal:", k_opt))

	return (k_opt)
}

approxNormal <- function(samples, samples_norm) {
	df_out <- data.frame(matrix(ncol=5, nrow=0,
							  dimnames=list(NULL, c("dist", "min", "max", "mu", "sigma"))))

	# assume Gaussian mixture model
	if (isMixture(samples)) {
		message("--- assuming Gaussian mixture model")
		tryCatch({
			k_opt <- estimateNumClusters(samples_norm)
			if (length(k_opt) > 0 && k_opt >= 2) {
				# compute actual values with optimal k
				message(paste("--- computing Gaussian mixture model with k =", k_opt))
				gmm <- normalmixEM(samples, k=k_opt, maxit=4048, maxrestarts=10, fast=TRUE)

				# add values to df
				data_min = min(samples)
				data_max = max(samples)
				for (k in 1:k_opt) {
					df_out[k,] <- c("norm", data_min, data_max, gmm$mu[k], gmm$sigma[k])
				}

				gmm.plot <- create_plot(gmm$x, gmm$mu, gmm$sigma, gmm$lambda)

				return (list(df_out, gmm.plot))
			}
		}, error = function(e) { message(paste("---- failed with error:", (e$message))) }
		)
	}

	message("--- assuming singular Gaussian model")
	# default to single Gaussian distribution
	mu <- mean(samples)
	sigma <- sd(samples)
	df_out[1,] <- c("norm", min(samples), max(samples), mu, sigma)

	nplot <- create_plot(samples, mu, sigma)

	return (list(df_out, nplot))
}

create_hist <- function(samples) {
	breaks <- pretty(range(samples), n = nclass.FD(samples), min.n = 1)
	bwidth <- breaks[2]-breaks[1]

	df <- data.frame(x = samples)
	plt <- ggplot(data=df) +
	  	   geom_histogram(aes(x, ..density..), binwidth = bwidth, colour = "black", fill = "white") +
		   xlab("Values") +
		   ylab("Density")

	return (plt)
}

create_plot <- function(samples, mu, sigma, lambda=1) {
	breaks <- pretty(range(samples), n = nclass.FD(samples), min.n = 1)
	bwidth <- breaks[2]-breaks[1]

	df <- data.frame(x = samples)
	plt <- ggplot(data=df) +
	       geom_histogram(aes(x, ..density..), binwidth = bwidth, colour = "black", fill = "white")

	for (k in 1:length(mu)) {
		plt <- plt + stat_function(geom = "line", fun = plot_mix_comps,
				  				   args = list(mu[k], sigma[k], lam = lambda[k]),
								   lwd = 1.5) 
	  }

	plt <- plt + xlab("Values") + ylab("Density")

	return (plt)
}

plot_mix_comps <- function(x, mu, sigma, lam) {
	  lam * dnorm(x, mu, sigma)
}

estimate_parameters <- function(file) {
	filename <- basename(file)
	output_dir <- paste0(filename, "_out")

	dir.create(file.path(getwd(), output_dir))

	message(paste("reading file", file))
	# TODO: add other file types
	df <- read.csv(file, stringsAsFactors=TRUE)

	out <- computeGMMs(df)

	df_out <- out[[1]]
	plots <- out[[2]]


	message(paste("writing output to", output_dir))
	write.csv(df_out, paste0(output_dir, "/", "distributions.csv"), row.names = FALSE)

	for (c in names(plots)) {
		ggsave(filename = paste0(output_dir, "/", c, ".pdf"),
			   plot = plots[[c]])
	}
}

# TODO: add for every file loop
# read spreadsheet
file <- "if/COVID-19_merged.csv"  # only for testing

