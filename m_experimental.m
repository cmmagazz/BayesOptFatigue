%% Experimental testing script.
% Simplified for a generic fatigue test.
%
% The aim is to run these section by section - pressing ctrl+shift+enter
% CMM 2020

resultsfilename='fatigueresults.mat';

%LOAD PREVIOUS DATA:
try
    load(resultsfilename)
catch
    disp('New dataset')
end
%% Sample Details:

sample='b'; %input sample name
try 
    exists(totalresults.(sample))
    disp('Continuing test')
catch
    disp('New sample')
end

%% Testing details
% you can break this section while waiting for it to complete. simply
% re-run with ctrl+enter to resume

testdet.protocol    = 'bayes step'; %choose a protocol from key: 
%                 stress step
%                 staircase
%                 probit
%                 life
%                 bayes staircase
%                 bayes step

%===== Necessary details
testdet.runout          = 1e6; %what's the runout value in cycles
testdet.step            = 50; %IF STEP TESTINGminimum stress step in MPa
testdet.startingstress  = 100; %minimum starting stress in MPa
totalresults.(sample).testdet=testdet;
testdet.prevsampdet     = f_findprevsamp(totalresults,testdet); %if there was a previous sample at that life (string or NaN)
testdet.testfreq        = 20400; %test frequency in Hz

%===== If using Bayesian methods, details of the prior space
testdet.theta=100:10:1000;
% set upper and lower bound of possible values of sigma here
testdet.sigma=1:10:100;

%===== any other details that are useful but not necessary
testdet.time        =   clock; %decimal time of experiment
testdet.torqwrench  =   4; %torque wrench setting in Nm

%put it all in the totalresults struct
totalresults.(sample).testdet=testdet;

totalresults.(sample)=f_testing(totalresults.(sample));

% End of testing:
close all
save(resultsfilename, 'totalresults')

%% SN Curve - all results
%First generate SN results for all the tests completed
totalresults=f_SNresults(totalresults);
p_SN(totalresults.SN.results)

%% EOS:
close all
save(resultsfilename, 'totalresults')


%% Prior plots
lprior=g_calcprior(totalresults.SN.results,testdet.theta,testdet.sigma,2);
theta=testdet.theta;
% set upper and lower bound of possible values of sigma here
sigma=testdet.sigma;
p_HPD(lprior)


