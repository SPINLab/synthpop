library('stringi')

# Author: Xander Wilcke
# Email: w.x.wilcke@vu.nl
#
# R script to regenerate dataset from distribution parameters. 
# Lorem Ipsum is used for natural language and other text values.

sampleNorm <- function(param) {
	rdata <- c()  # empty vector

	# find column indices
	cnames <- colnames(param)
	mu_idx <- which(cnames=="mu")
	sigma_idx <- which(cnames=="sigma")

	num_distinct <- rep(floor(unique(param$num_distinct)/nrow(param)),nrow(param))  # per distribution
	mod_distict <- unique(param$num_distinct)%%nrow(param)
	i <- 1
	while (mod_distict > 0) {
		num_distinct[i] <- num_distinct[i] + 1
		mod_distict <- mod_distict - 1

		i <- i+1
	}

	# ugly fix to remedy impossible situation in which num_distinct < nrow(param)
	# this should be fixed in parameter estimation script
	num_distinct <- num_distinct[num_distinct > 0]

	len_distinct <- length(num_distinct)

	num_values <- rep(floor(unique(param$num_values)/len_distinct),len_distinct)  # per distribution
	mod_values <- unique(param$num_values)%%len_distinct
	i <- 1
	while (mod_values > 0) {
		num_values[i] <- num_values[i] + 1
		mod_values <- mod_values - 1
		
		i <- i+1
	}

	for (i in 1:len_distinct) {
		cdata <- rnorm(n=num_values[i],
		               mean=param[i,mu_idx],
			           sd=param[i,sigma_idx])

		if (num_distinct[i] == 1){
			cdata_sample <- sample(cdata, size=1)
			rdata <- c(rdata, rep(cdata_sample, num_values[i]))
		} else if (num_distinct[i] < num_values[i]) {
			cdata_sample <- sort(sample(cdata, size=num_distinct[i], replace=FALSE))

			psample <- c()
			for (j in 1:num_distinct[i]) {
				if (j == 1) {
					psample <- c(psample, pnorm(cdata_sample[j], param[i,mu_idx], param[i,sigma_idx]))
				} else {
					pb <- pnorm(cdata_sample[j-1], param[i,mu_idx], param[i,sigma_idx])
					pe <- pnorm(cdata_sample[j], param[i,mu_idx], param[i,sigma_idx])
					psample <- c(psample, (pe-pb))
				}				
			}

			psample <- psample * 1/sum(psample)  # scale with sum 1
			rdata <- c(rdata, sample(cdata_sample, size=num_values[i], replace=TRUE, prob=psample))
		} else {
			rdata <- c(rdata, cdata)
		}
	}

	rdata <- sample(rdata)

	return (rdata)
}

sampleUnif <- function(param) {
	num_values <- unique(param$num_values)
	num_distinct <- unique(param$num_distinct)

	if (num_distinct == 1) {
		rdata <- rep(runif(1, min=param$min, max=param$max), num_values)
	} else if (num_distinct < num_values) {
		rdata_diff <- param$max - param$min
		rdata_incr <- rdata_diff / (num_distinct - 1)
		rdata_values <- c(param$min)

		v <- param$min
		for (i in 1:(num_distinct - 1)) {
			v_new <- v + rdata_incr
			rdata_values <- c(rdata_values, v_new)

			v <- v_new
		}

		p <- rep(1.0 / num_distinct, num_distinct)  # add some difference in counts
		rdata <- sample(rdata_values, size=num_values, replace=TRUE, prob=p)
	} else {
		rdata <- runif(param$num_values, min=param$min, max=param$max)
	}

	return (rdata)
}

generate_factors <- function(param, cdata) {
	terms_min <- unique(param$terms_min)
	terms_max <- unique(param$terms_max)

	if (is.na(terms_min) || is.na(terms_max)) {
		return (cdata)
	}

	keys <- unique(cdata)
	key_sentence_map <- sapply(keys, int_to_sentence,
							   terms_min=terms_min, terms_max=terms_max) 

	for (key in keys) {
		cdata[which(cdata==key)] <- key_sentence_map[as.character(key)]
	}

	return (cdata)
}

int_to_sentence <- function(key, vec=c(), terms_min, terms_max) {
	vec[as.character(key)] <- generate_sentence(terms_min, terms_max)

	return (vec)
}

generate_sentence <- function(terms_min, terms_max) {
	num_terms <- sample.int(terms_max-terms_min+1, size=1) + terms_min - 1

	terms <- c()
	while (length(terms) < num_terms) {
		lipsum <- stri_rand_lipsum(nparagraphs=1, start_lipsum=FALSE)
		terms <- c(terms, strsplit(lipsum, ' ')[[1]])
	}

	sentence <- paste(terms[1:num_terms], collapse=' ')
	if (terms_max <= 1) {
		# assume identifier
		sentence <- toupper(sentence)
	} else {
		# assume natural language
		end <- substring(sentence, nchar(sentence))
		if (end != '.') {
			if (end == ',') {
				substring(sentence, nchar(sentence)) <- '.'
			} else {
				sentence <- paste0(sentence, '.')
			}
		}
	}

	return (sentence)
}

generate_column <- function(param, lvl, num_values) {
	df <- data.frame()

	# sample from distributions
	# ignore exceeding min/max values are we're dealing with estimates anyway
	distribution <- unique(param$dist)
	if (distribution=="gmm") {
		cdata <- sampleNorm(param)
	} else {
		cdata <- sampleUnif(param)
	}
	
	# cast to integers if type is non numeric
	dtype <- unique(param$type)
	if (!(dtype=="numeric")) {
		distinct <- unique(cdata)
		for (i in 1:length(distinct)) {
			cdata[cdata==distinct[i]] <- i
		}
		cdata <- as.integer(cdata)
	}

	# convert to factor data
	if (dtype=="factor") {
		cdata <- generate_factors(param, cdata)
	}

	# add NAs if needed
	if (length(cdata) < num_values) {
		filler <- rep(NA, num_values-length(cdata))
		cdata <- sample(c(cdata, filler))  # merge and shuffle
	}

	if (length(df) <= 0) {
		# initialize dataframe
		df <- data.frame(cdata)
		colnames(df) <- lvl
	} else {
		df[lvl] <- cdata
	}

	return (df)
}

generate_dataset <- function(filename, columnwise=FALSE) {
	output_dir <- dirname(filename)
	message(paste("output directory set to", output_dir))

	message(paste("reading file", filename))
	params <- read.csv(filename, stringsAsFactors=TRUE)

	df <- data.frame()
	num_values <- max(params$num_values)
	for (lvl in levels(params$name)) {
		message(paste(" generating column", lvl))
		param <- params[params$name==lvl,]

		df_out <- generate_column(param, lvl, num_values)

		if (columnwise) {
			name_out <- paste0(output_dir, "/", "dataset_", lvl, ".csv")
			write.csv(df_out, name_out, row.names = FALSE)
		} else {
			if (length(df) <= 0) {
				# initialize dataframe
				df <- data.frame(df_out)
				colnames(df) <- lvl
			} else {
				df[lvl] <- df_out
			}
		}
	} 

	if (!columnwise) {
		name_out <- paste0(output_dir, "/", "dataset.csv")
		write.csv(df, name_out, row.names = FALSE)
	} else {
		message("produce a single file with `$ paste -d ',' dataset_*.csv > dataset.csv`")
	}
}

#csv_file <- "COVID-19_merged.csv_out/distributions.csv"

