function y= g_UTILITY(stress,theta,sigma,lprior)
% Utility function for Gearhart method
% Conceptually: we want to find a stress where, whether it breaks or fails,
% it will give us the highest information gain. This is also scaled by how
% likely it is, given our previous data, for us to collect this datapoint. 
%
% input:
%   stress: a double to compute the probabilities at
%   theta: space of theta values allowed, vector
%   sigma: space of sigma values allowed, vector
%   lprior: prior distribution across the theta/sigma space
%
% output:
%   y: value for utility at this stress, given this prior distribution
%
% CMM, RR, CG 2021

%Calculate the joint posterior, and the probability of breaking or not 
% breaking at this stress. g_calcprior can be used to calculate this with
% "fake" data. 

[lpostbreak,xbreak]=g_calcprior([stress,NaN,NaN],theta,sigma,lprior);
[lpostnotbreak,~,~,xnotbreak]=g_calcprior([NaN,NaN,stress],theta,sigma,lprior);
xnotbreak=1-xnotbreak; %x2surv is the cdf of failing at runoutstress, so 
%need to flip to make it survive

%calculate the marginal probability
probbreak=exp(lprior).*xbreak;
probnotbreak=exp(lprior).*xnotbreak;
%integrate
probbreak=sum(sum(probbreak));
probnotbreak=sum(sum(probnotbreak));

%normalise the posterior, and finally the integral across the prior:
normbreak=sum(sum(exp(lpostbreak)));
lpostbreak=lpostbreak-ones(size(lprior)).*log(normbreak);
normnotbreak=sum(sum(exp(lpostnotbreak)));
lpostnotbreak=lpostnotbreak-ones(size(lprior)).*log(normnotbreak);
%information is p*logp
infobreak=sum(sum(lpostbreak.*exp(lpostbreak)));
infonotbreak=sum(sum(lpostnotbreak.*exp(lpostnotbreak)));
%output: negative for original convenience of finding the minimum of y
%(which is the maximum of utility)
y=-(probbreak*infobreak+probnotbreak*infonotbreak);

end

