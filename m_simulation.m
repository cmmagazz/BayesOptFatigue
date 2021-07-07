%% Simulation Main Script
% This script is designed to run in sections. 
% The first section initialises the variables, organised into two primary
% structs: TestingSet and ResultSet. 
% TestingSet contains the simulated sample dataset
% ResultSet contains the results of testing this simulated dataset

%=============================================
clear all
close all
%=============================================
%Testing Set parameters and setup
TestingSet.details.basquin.c=350;
TestingSet.details.basquin.alpha=500;
TestingSet.details.basquin.beta=-log(1/6)/log(6);
TestingSet.details.N=linspace(0.5,10,1000); % in log N
TestingSet.details.DAMq=1;%damage model: 1=none, 2=miner's rule
%Sample distribution types - a 2 element vector with:
%       first  element: sample distribution: 1=gaussian, 2=weibull
%       second element: mean model: 1=basquin curve, 2=bi-linear curve
TestingSet.details.SDistq=[1,2]; %sample distribution: 1=gaussian, 2=weibull

TestingSet.details.numsamp=20; %number of specimens
TestingSet.details.width=10; %width in sigma

TestingSet.meanFS=f_createsample(TestingSet.details);
%=============================================
%ResultSet parameters
ResultSet.details.startingstress=150;
ResultSet.details.runout=6; %set the runout value
ResultSet.details.step.stepsize=100;

%=============================================
%constant life prior constants 
maxtheta=700;
mintheta=300;
theta=mintheta:2:maxtheta;
% set upper and lower bound of possible values of sigma here
maxsigma=100;
minsigma=1;
sigma=minsigma:2:maxsigma;

ResultSet.details.theta=theta;
ResultSet.details.sigma=sigma;
%% Simple run
close all
ResultSet.details.protocol = 'bayes staircase';
ResultSet.raw=f_testingProtocol(TestingSet,ResultSet);

%% Plotting scripts: to be run line by line if desired
figure()
plot(ResultSet.raw.shannon)
title('Evolution of Shannon Information')
xlabel('Sample Number')
ylabel('Shannon Information')

%Plot the staircase plot of results
p_staircaseplot(ResultSet.raw.failurestress)


%Plot the HPD
p_HPD(ResultSet.raw.lprior)


