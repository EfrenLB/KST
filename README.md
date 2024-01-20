This repository contains the script for quantifying maximum differences in cumulative distribution functions (D- and P-values from Kolmogorov–Smirnov tests) between different sample populations (the pan-Arctic domain and INTERACT stations, with or without Russian stations) across eight CMIP6 ESMs and eight ecosystem variables. The code randomly samples the same number of grid cells from all ESMs, equal to the smallest population size among all models (CanESM5 with 496 data points, excluding ocean and Greenland Ice Sheet pixels). To minimize potential artifacts from the sample size choice, 100 replicates of random sample populations of 496 data points per ESM and variable are retrieved. Additionally, the script also extracts the first (25%), second (median), and third (75%) quartile (Q1-3) values of the distribution functions for each ESM and ecosystem variable for the same populations. 

To execute the script, use the KS-test.R script to calculate D- and P-values, and Q1-3. Additional details about the rationale behind the script and its purpose are available in the scientific paper: López-Blanco, E., Topp-Jørgensen, E., Christensen, T. R., Rasch, M., Skov, M., Arndal, M. F., Bret-Harte, S. M., Callaghan, T. V., and Schmidt, N. M.  "Towards an increasingly biased view on Arctic change." (2024) Nature Climate Change, https://www.nature.com/articles/s41558-023-01903-1.


