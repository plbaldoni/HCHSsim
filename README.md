# HCHSsim
Set of R codes for the simulation study using HCHS survey sampling design

## Overview

This set of R codes replicate the results from the simulation study implementing the bootstrap and multiple imputation correction of standard errors from the health outcome model with biomarker calibrated nutrients. The codes presented here generate a target population using a three-stage design, draw samples from it using a stratified complex survey sampling scheme, and implement the aforementioned methods for standard error correction.

## Required software

The entire pipeline is written in R. All necessary packages are loaded in the headers of the codes.

## How-to

To run the simulation pipeline, you should:

1. Run *generateTargetPopulation.R*. This will generate a list of participants using a three-stage design with stratification, block groups, and households.
2. Run *generateTargetPopulationData.R*, populate the target population data sets with several variables to be used in this simulation study.
3. Run *sampleTargetPopulation.R*, which will generate 1000 samples from the target population using a stratified three-stage sampling scheme. Sampling probability weights and design variables are included in the final data sets.
4. Run *runBootstrap.R* and *runMI.R* to run the proposed approaches for standard error correction.

The *runRaking.R* script attempts to implement the raking estimator as a third alternative to correct the standard errors from the outcome model with biomarker-calibrated nutrients as the main exposure. The current implementation needs to be revised.
