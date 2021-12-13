library('stringi')

sampleNorm <- function(param) {
	rdata <- c()  # empty vector

	# find column indices
	cnames <- colnames(param)
	mu_idx <- which(cnames=="mu")
	sigma_idx <- which(cnames=="sigma")

	num_values <- unique(param$num_values)/nrow(param)  # per distribution
	for (i in 1:nrow(param)) {
		rdata <- c(rdata, rnorm(n=floor(num_values),
							    mean=param[i,mu_idx],
							    sd=param[i,sigma_idx]))
	}

	# compensate for rounding errors on num_values
	if (length(rdata) < unique(param$num_values)) {
		j = sample.int(nrow(param), size=1)
		rdata <- c(rdata, rnorm(n=1,
								mean=param[j,mu_idx],
								sd=param[j,sigma_idx]))
	}

	# shuffle data to avoid ordered distributions from loop
	rdata <- sample(rdata)

	return (rdata)
}

sampleUnif <- function(param) {
	rdata <- runif(param$num_values, min=param$min, max=param$max)

	return (rdata)
}

resample <- function(param, rdata) {
	num_values <- unique(param$num_values)
	num_distinct <- unique(param$num_distinct)

	# resample from a distinct subset
	if (num_distinct < num_values) {
		rdata_sample = sample(rdata, size=num_distinct, replace=FALSE)
		rdata <- sample(rdata_sample, size=num_values, replace=TRUE)
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
		cdata <- sapply(cdata, round)  # use round to avoid hard cut
		cdata <- as.integer(cdata)
	}

	# resample data to match the number of distinct values
	cdata <- resample(param, cdata)

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

