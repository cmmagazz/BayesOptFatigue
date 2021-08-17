function [beststress,lprior,shannon]=B_simulate(failurestress,theta,lpriorq,lprior)
%Take the failuretally of the samples so far (if there are any) and find 
%the next stress to test. 
%
%INPUTS:
%   failurestress: the failure data of the current results
%   theta: the sample points for the mean failure stress (e.g
%       theta=mintheta:0.5:maxtheta)
%   sigma: the sample points for the std of the failure stress (e.g
%       sigma=minsigma:0.5:maxsigma)
%OUTPUTS:
%       beststress: the optimal next stress to test
%       lprior: the log of the prior, as an array
%       shannon: the shannon information (=p*log(p))

[lprior,~,shannon]=g_calcprior(failurestress,theta,lprior,lpriorq);
beststress=g_bayesbeststress(theta,lprior);
end