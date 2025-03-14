# Temporal M-quantile models and robust bias-corrected small area predictors

### Bugallo M. [1], Morales D. [1], Salvati N. [2] and Schirripa Spagnolo F. [2]

### Center of Operations Research, Miguel Hernández University of Elche, Spain [1]
### Department of Economics and Management, University of Pisa, Italy [2]

## Objective


This repository contains R scripts for implementing the new statistical methodology based on Time-Weighted M-Quantile (TWMQ) linear models in a practical case. Simulated data is included to illustrate its application and analysis. The corresponding manuscript is published in the open-access archive *arXiv* and can be accessed at: [http://arxiv.org/abs/2407.09062](http://arxiv.org/abs/2407.09062).

## Installation and requirements

The R code available in this repository was tested on a 3.3 GHz Intel Xeon W Mac computer with 32 GB RAM. For data manipulation and task processing, the following R packages are required:

-   library(MASS)
-   library(fpp2) 
-   library(dplyr)
-   library(fastDummies)
-   library(agricolae)


To facilitate code manipulation, we recommend using the **Graphical User Interface (GUI)** [RStudio](https://posit.co/downloads/). Running the scripts and reproducing the results is straightforward. Once all the required libraries are installed, simply download the project sources and execute the main script ***TWMQmodelsUsageExample.R***.

## Contents

1. **Main script – *TWMQmodelsUsageExample.R***: Implements the TWMQ models with both continuous and categorical covariates, computes the predictors, and estimates their mean squared error (MSE).  

2. **Simulated data**: Test dataset for applying the methodology, simulated in the script ***TWMQmodelsUsageExample.R***.

3. **Auxiliary scripts – *QRLM.R* and *QRLMweights.R***: Contain essential functions for model fitting, temporal predictor computation, and the proposed MSE estimations.  
