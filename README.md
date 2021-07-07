# BayesOptFatigue

## Overview

The repo contains a set of functions for performing bayesian statistics on fatigue data where for each sample the stress tested and whether the sample failed by a given number of cycles. The Functions can take either staircase-type (one stress per sample) or step-tpe (multiple stresses per sample increasing stress by a step size until failure).
This is a companion to the paper Bayesian Optimised Collection Stretegies for fatigue Testing: Constant Life Testing (https://arxiv.org/abs/2107.02685)

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

![Flow Chart of Code](./code_flow_chart.png)

## Functions

Main Scripts:

* **m_simulation** - running protocols on simulated datasets
* **m_experimentaltesting** - running protocols on real datasets during testing

Calculation Functions:

* **f_testing_protocol** - manages protocols based on previous samples
* **f_createsample** - creates simulated sample set
* **f_testingDAM** - tests simulated samples based on their determenistic SN curve
* **B_simulate** - runs constant life bayesian staircase protocol
* **B_STEP_simulate** - runs constant life bayesian step protocol
* **g_bayes_beststepsize_stepstart** - runs step utility function for all possible test parameters and optuts maximum value
* **g_bayesbeststress** - runs staircase utility function for all possible test parameters and optuts maximum value
* **g_calcprior** - updates parameter probabilities based on data and bayes theorem and outputs new log(prior)
* **g_STEP_UTILITY** - calculates the step utility function for a set of protocol parameters and prior model distribution
* **g_UTILITY** - calculates the staircase utility function for a set of protocol parameters and prior model distribution
* **f_HPD** - calculates the heighest prosterior density interval of a joint posterior distribution

Plotting Functions:

* **p_HPD** - Plots the joint posterior HPD interval
* **p_staircaseplot** - Plots staircase plots for data

## Key Variables

* **TestingSet** - Details of testing
* **ResultsSet** - Results of testing so far
* **Theta** - Vector of allowed means
* **Sigma** - Vector of allowed standard deviations
* **lprior** - joint posterior/prior of data, logged for storage
* **failurestress** - Results of tests - Failure stress (MPa), Numebr of steps to failure, Runout Stress (MPa), Cyces for runout life (log(cycles)), Failure Tally (1=fail, 0=runout)
* **stepsize** - step size used for step or staircase testing
* **startingstress** - starting stress for all protocols and test stress for life testing
* **Shannon** - shanon information of joint posterior





