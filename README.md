# BayesOptFatigue

## Overview

The repo contains a set of functions for performing bayesian statistics on fatigue data where for each sample the stress tested and whether the sample failed by a given number of cycles. The Functions can take either staircase-type (one stress per sample) or step-tpe (multiple stresses per sample increasing stress by a fixed step size until failure).
This is a companion to the paper Bayesian Optimised Collection Strategies for fatigue Testing: Constant Life Testing (https://arxiv.org/abs/2107.02685)

Based on the data, the functions can:

* Plot the results on standard plots
* Perform Bayesian analysis to model the data with Normal, Log-Normal, Weibull distributions
* Compare evidence for each model using Bayes Information Criterion
* Perform Bayesian Maximum Entropy sampling to determine the optimal next test location and protocol parameters

Functions are also supplied for creating simulated sample sets based on Monte-Carlo sampling with a given input distribution.

Supported Protocols are:

* Life Testing
* Staircase
* Stress Step
* Probit
* Bayesian Staircase
* Bayesian Step

A schematic flow diagram of the software is shown below:

![Flow Chart of Code](./code_flowchart.png) 

## Functions

Main Scripts:

* **m_simulation** - running protocols on simulated datasets
* **p_BIC_convergence** - runs many protocols on simulated data in paralell to demonstrate effects on distinguishability between models which changing variables and plots results.
* **m_experimental** - running protocols on real datasets during testing

General Calculation Functions:

* **Simulation framework**
* f_testing_protocol - tests a simulated sample using a set protocol
* f_createsample - creates simulated sample set
* f_testingDAM - tests simulated samples based on their determenistic SN curve
* **Experimental testing framework**
* f_testing - tests a sample using a set protocol
* f_findprevsamp - read through the experimental campaign variable to find previous completementary samples
* f_SNresults - compile the results of an experimental campaign

Bayesian Calculation Functions
* **g_bayes_beststepsize_stepstart** - runs step utility function for all possible test parameters and optuts maximum value
* **g_bayesbeststress** - runs staircase utility function for all possible test parameters and optuts maximum value
* **g_calcprior** - updates parameter probabilities based on data and bayes theorem and outputs new log(prior)
* **g_STEP_UTILITY** - calculates the step utility function for a set of protocol parameters and prior model distribution
* **g_UTILITY** - calculates the staircase utility function for a set of protocol parameters and prior model distribution
* **f_HPD** - calculates the heighest prosterior density interval of a joint posterior distribution
* **g_calcaic** - calculates AIC, BIC, and log-likelihood functions for normal, log-normal, and weibull distributions
* **g_calcloglike** - calculates liklelihood functions and macimum log-likelihood

Plotting Functions:

* **p_HPD** - Plots the joint posterior HPD interval
* **p_staircaseplot** - Plots staircase plots for data
* **p_contourHPD** - Plots joint posterior HPD as a contour plot
* **p_priorcomparison** - Plots comparison between multiple joint posteriors and includes marginal parameter probabilities
* **p_SN** - Plots an SN curve with (if available) three points: penultimate runout, final failure stress, and cycles to failure
* **c_convergence** - Plots convergence, or error comparison, diagrams as in the paper linked. 


## Key Variables

* **TestingSet** - Details of testing
* **ResultsSet** - Results of testing so far
* **Theta** - Vector of allowed means and scales
* **lprior** - joint posterior/prior of data, logged for storage
* **failurestress** - Results of tests - Failure stress (MPa), Numebr of steps to failure, Runout Stress (MPa), Cyces for runout life (log(cycles)), Failure Tally (1=fail, 0=runout)
* **stepsize** - step size used for step or staircase testing
* **startingstress** - starting stress for all protocols and test stress for life testing
* **Shannon** - shanon information of joint posterior


## Dependencies
* ShadedErrorBar from https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
* MATLAB Parallel Computing Toolbox
* MATLAB Statistics and machine learning toolbox
* MATLAB Global optimisation toolbox (suggested)

## Requirements
* MATLAB version 2019a or later



