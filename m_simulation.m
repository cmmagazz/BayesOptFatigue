%% Simulation Main Script
% This script is designed to run in sections. 
% The first section initialises the variables, organised into two primary
% structs: TestingSet and ResultSet. 
% TestingSet contains the simulated sample dataset
% ResultSet contains the results of testing this simulated dataset

%=============================================
clear all
close all
clc
%=============================================
%Testing Set parameters and setup
TestingSet.details.basquin.c=350;
TestingSet.details.basquin.alpha=500;
TestingSet.details.basquin.beta=-log(1/6)/log(6);
TestingSet.details.N=linspace(0.5,10,1000); % in log N
TestingSet.details.DAMq=1;%damage model: 1=none, 2=miner's rule

%Sample distribution types - a 2 element vector with:
%       first  element: sample distribution: 1=gaussian, 2=weibull
%       optional second: mean model: 1=basquin (default), 2=bi-linear
TestingSet.details.SDistq=[1,2]; 

TestingSet.details.numsamp=100; %number of specimens
TestingSet.details.width=20; %width in sigma

TestingSet.meanFS=f_createsample(TestingSet.details); %create the monte-carlo sample set
%=============================================
%Constant life prior constants: determine the space of parameters: 
% Insert range of parameters, see function for details
% Choose from: 'norm','lognorm','2pwbl','3pwbl'
% Choose numel 
% ResultSet.details=f_setupresultsdist([[300, 700];[1, 150]],'norm',150);
% ResultSet.details=f_setupresultsdist([[300, 700];[1, 150]./150],'lognorm',150);
% ResultSet.details=f_setupresultsdist([[300, 700];[1, 150]],'2pwbl',150);
ResultSet.details=f_setupresultsdist([[300, 700];[1, 150];[0,300]],'3pwbl',50);

%=============================================
%ResultSet parameters
ResultSet.details.startingstress=150; %starting stress in MPa 
ResultSet.details.runout=6; %set the runout value in log(cycles)
ResultSet.details.step.stepsize=100; %step size in MPa
%% Simple run
%Select a protocol from
%     'stress step'
%     'staircase'
%     'probit'
%     'life'
%     'bayes staircase'
%     'bayes step'
clc
close all
ResultSet.details.protocol = 'bayes step';
%Plot a diagnostic SN curve during collection? 1=yes, 0=no
ResultSet.plotq=0;
%Run the simulated test
ResultSet=f_testingProtocol(TestingSet,ResultSet);

%% Plotting scripts: to be run line by line if desired
%Shannon information only available for bayesian protocols ** see below
figure()
plot(ResultSet.raw.shannon)
title('Evolution of Shannon Information')
xlabel('Sample Number')
ylabel('Shannon Information')

%% Plot the staircase plot of results
p_staircaseplot(ResultSet.raw.failurestress)


%% Plot the HPD
p_HPD(ResultSet.raw.lprior)

%% Evaluate the joint posterior from a dataset
% ** see here to evaluate a 
[lprior,~,shannon]=g_calcprior(ResultSet.raw.failurestress,ResultSet.details.theta,ResultSet.details.sigma);
p_HPD(lprior)

%% Run through 3d prior
for i=1:size(ResultSet.raw.lprior,3)
    pcolor(exp(squeeze(ResultSet.raw.lprior(:,:,i))));
    pause(0.1)
    colorbar
    hold on
end