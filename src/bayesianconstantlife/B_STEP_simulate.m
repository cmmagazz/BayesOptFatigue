function [beststress,beststep,lprior,shannon]=B_STEP_simulate(failurestress,theta,startingstress,minstressstep,varargin)
%Take the failuretally of the samples so far (if there are any) and find 
%the next stress to test. 
%
%INPUTS:
%   failurestress: vector with layout [failure stress,~,runout stress]
%   thetaR: the sample points for the mean failure stress (e.g
%       theta=mintheta:0.5:maxtheta)
%   sigmaR: the sample points for the std of the failure stress (e.g
%       sigma=minsigma:0.5:maxsigma)
%OUTPUTS:
%   beststress: the optimal next starting stress to test
%   beststep: the optimal next stress step to test
%Calculate your prior
if numel(varargin)==1
    [lprior,~,shannon]=g_calcprior(failurestress,theta,varargin{1});
elseif numel(varargin)==0
    [lprior,~,shannon]=g_calcprior(failurestress,theta);
end
[beststress,beststep]=g_bayes_beststepsize_stepstart(theta,lprior,minstressstep,startingstress);
end