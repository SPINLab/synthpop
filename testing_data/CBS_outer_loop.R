library(haven)
library(labelled)

setwd("H:/accessible_file_statistics")

#log <- file("logmain.txt", open="a")
#sink(log, type="message")
source("estimateParameters.R")

files <- file("accessible_file.txt", 'r')
while ( TRUE ) {
	f = readLines(files, n=1)
	if ( length(f) == 0 ) {
		break
	}

	# skip comments
	if ( substring(f, 1, 1) == "#" ) {
		next
	}
	
	# rmv drive letter and split remainder per directory
	f_no_drive = substring(f, 4, nchar(f))
	split = strsplit(f_no_drive, "\\\\")[[1]]
	
	# basename
	f_name = split[length(split)]
	
	if ( !file.exists(f)) {
		print(paste("File not found:", f_name))
		next
	}
	
	print(paste("Processing file:", f_name))	
	ext = tolower(substring(f_name, nchar(f_name)-2, nchar(f_name)))
	if ( ext == "sav" ) {
		df = as.data.frame(read_sav(f), stringsAsFactors=TRUE)
	} else {
		print("Skipping file")
		next
	}

	tryCatch({
		# replace haven labels by actual values
		df <- labelled::unlabelled(df)

		# check which columns are numeric and set type
		for (i in 1:length(colnames(df))) {
			if (all(grepl("^[-,\\.0-9]+$", df[,i]))) {
				df[,i] <- as.numeric(df[,i])
				class(df[,i]) <- "numeric"
			}
		}

		estimate_params(df, paste(f_name, "_out"))
	}, error = function(e) { message(paste("----- failed with error:", (e$message))) }
	)

	rm(df)
	gc()
}
close(files)
print("Done")
