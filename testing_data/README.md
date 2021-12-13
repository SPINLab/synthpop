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

1) Start R (or RStudio) and source the parameters estimation script:

`> source("estimate_parameters.R")`

2) Call the parameter estimation function with a file descriptor:

`> estimate_parameters("/path/to/file")`

3) Once the above function has completed (this can take a few hours), source the data generation script:

`> source("generate_mock_data.R")`

4) Call the data generation function with a file descriptor to the just-estimated parameters:

`> generate_dataset("/path/to/parameters.csv")`

