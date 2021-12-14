# Synthetic Data Generation by Statistical Inferencing 

Synthetic data is becoming an increasingly important topic within and outside of academic research, especially when privacy sensitive data is being considered. Directly using privacy sensitive data in experiments or simulations is generally considered unsafe, as there is a non-zero chance that some elements of these data find their way into the results. To generate such a synthetic dataset it is often necessary to have access to the original and privacy-sensitive data. This might not always be possible, however, for example when the providers of the data have restricted access to the dataset or when this dataset is only accessible through a restricted environment. To reduce the dependency on these data, it is sensible to do most development and testing on a mock dataset which looks and feels like the original dataset, but of which the content is entirely fabricated.

The scripts published here can be used to generate mock data that have a distribution and range that is similar to the original data. This is achieved using a two-step approach: first, the parameters of the underlying distributions are estimated column by column from the original data, thereby combining statistical inferencing and Monte Carlo sampling, after which these parameters are used to locally reconstruct an approximation of these data. See the PDF for a detailed overview of the method.

## Prerequisites

- R or RStudio
- The following R libraries:
	- fitdistrplus
	- FNN
	- ggplot2
	- mixtools
	- stringi

## HowTo

1) Read the original data into R (of RStudio) via `read.csv()` or similar:

`> data <- read.csv("path/to/file.csv")`

2) Source the parameters estimation script and call the `estimate_parameter()` function with, optionally, the destination directory (default `output/`) and the number of samples to estimate from (default `10e4`):

`> source("estimate_parameters.R")`

`> estimate_parameters(data[, <output dir>][, <number of samples>])`

This can take some time to finish, depending on the number of samples and the size of the data. On completion, the output directory contains the file `distributions.csv`, in which each row holds the parameters of a single distribution. How well these distributions fit the data can be seen in the accompanied plots.

3) Generate mock data by sourcing the relevant script and by calling the `generate_dataset()` function. Set `columnwise=TRUE` if keeping the whole dataset in memory is unfeasible:

`> source("generate_mock_data.R")`

`> generate_dataset("/path/to/distributions.csv"[, columnwise=TRUE/FALSE])`

