# Project Name

Supplementary material and code for "Uncertain beginnings: A comparison of the accuracy and precision of methods estimating extreme chornological events, by Alfredo Cortell-Nicolau, Enrico R. Crema and Alastair Key. Currently under review.

## Contents and structure
The present project contains the following four folders:

* Figures: Figures used in the manuscript and supplementary information.
* R: Scripts necessary for the full reproducibility of the paper (see 'usage' section).
* Results: Results after simulations, in rds format.
* Simu\_data: Simulated data used for the simulations. 

## Usage
To fully reproduce the article, the different scripts must be run in the order shown in the folder R, and below. Paths are local, so they should work out of the box, but please have a look if you need to modify them in your setup. Here's a brief summary of what each script does.

1. 01\_Functions.R: Generates custom functions, which are necessare both for the implementation of OLE and BPM through R-based Oxcal.
2. 02\_Simu\_dates: Generates the dates according to the different population models, accounting for the back-calibration process.
3. 03\_Combined\_sims.R: This is the main script for the paper. It runs all the different implementations for all the models and parameter combinations, on the dates produced with the script above. This script took several days to run parallellised at 35 cores (almost two weeks). We'd advice to try first with a small number of simulations before sending it to a cluster. We discourage its use if there's no access to parallel computing.
4. 04\_Data\_output.R: This is just a utility script that homogeneises some of the data outputs produce above, for their systematic use in the scripts below.
5. 05\_Results.R: Produces the plots and computations of the results for the main manuscript.
6. 06\_Results\_SI.R: Produces the plost and computations of the results for the supplementary information.

## Wrap up
And I guess that's all you need to know, but please do reach out if you have any doubt on how to implement this!




