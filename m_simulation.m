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
addpath(genpath('src'));
%=============================================
%Testing Set parameters and setup
TestingSet.details.basquin.c=300;
TestingSet.details.basquin.alpha=600;
TestingSet.details.basquin.beta=-log(1/6)/log(6);
TestingSet.details.N=linspace(0.5,10,1000); % in log N
TestingSet.details.DAMq=1;%damage model: 1=none, 2=miner's rule

%Sample distribution types - a 2 element vector with:
%       first  element: sample distribution: 1=gaussian, 2=weibull
%       optional second: mean model: 1=modified (default), 2=basquin,
%       3=bi-linear,4=wohler
TestingSet.details.SDistq=[1,1]; 

TestingSet.details.numsamp=1; %number of specimens
TestingSet.details.width=0; %width in sigma

TestingSet.meanFS=f_createsample(TestingSet.details,1); %create the monte-carlo sample set, 1=plot, 0=no plot

%=============================================
%Constant life prior constants: determine the space of parameters: 
% Insert range of parameters, see function for details
% Choose from: 'norm','lognorm','2pwbl','3pwbl','gumbell'
% Choose numel 
ResultSet.details=f_setupresultsdist([[200, 600];[1, 400]],'norm',[250,200]);
% ResultSet.details=f_setupresultsdist([[300, 700];[1, 150]./150],'lognorm',150);
% ResultSet.details=f_setupresultsdist([[300, 700];[1, 150]],'2pwbl',150);
% ResultSet.details=f_setupresultsdist([[300, 700];[1, 150];[0,300]],'3pwbl',50);

%=============================================
%ResultSet parameters

ResultSet.details.step.stepsize=20; %step size in MPa
ResultSet.details.startingstress=350; %starting stress in MPa 
ResultSet.details.runout=6; %set the runout value in log(cycles)

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
TestingSet.details.numsamp=60; %number of specimens

ResultSet.details.protocol = 'bayes staircase';
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
p_staircaseplot(ResultSet.raw.failurestress,'Staircase Plot',0)


%% Plot the HPD
ResultSet.raw.lprior=g_calcprior(ResultSet);
p_priorcomparison({ResultSet.raw.lprior},ResultSet,{'60 samples'},[0,100],[200,600])
%% Evaluate the joint posterior from a dataset
% ** see here to evaluate a 
[lprior,~,shannon]=g_calcprior(ResultSet);
p_HPD(lprior)

%% Calculate BIC
close all
[aic,bic]=g_calcaic(ResultSet.raw,1);
%% Run through 3d prior
for i=1:size(lprior,3)
    pcolor(exp(squeeze(lprior(:,:,i))));
    pause(0.1)
    colorbar
    hold on
end
%%
p_HPD(squeeze(sum(lprior,3)))
%% Another type of plot looking at prior over different numbers of samples
figure()
ResultSetT=ResultSet;
itotry=[10,30,60];
l=[];
for idx=1:length(itotry)
    ResultSetT.raw.failurestress=ResultSet.raw.failurestress(1:itotry(idx),:);
    l{idx}=g_calcprior(ResultSetT);
end
p_priorcomparison(l,ResultSet,{'10 Samples','30 Samples','60 Samples'},[1,150],[250,550])

%}