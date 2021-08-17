function y= g_UTILITY(stress,ResultSet)
% Utility function for Gearhart method
% Conceptually: we want to find a stress where, whether it breaks or fails,
% it will give us the highest information gain. This is also scaled by how
% likely it is, given our previous data, for us to collect this datapoint. 
%
% input:
%   stress: a double to compute the probabilities at
%   theta: space of parameter values allowed, cell array of vector
%   lprior: prior distribution across the theta/sigma space
%
% output:
%   y: value for utility at this stress, given this prior distribution
%
% CMM, RR, CG 2021

%Calculate the joint posterior, and the probability of breaking or not 
% breaking at this stress. g_calcprior can be used to calculate this with
% "fake" data. 
lprior=ResultSet.raw.lprior;

ResultSet.raw.failurestress=[stress,NaN,NaN];
[lpostbreak,xbreak]=g_calcprior(ResultSet,[],lprior);
ResultSet.raw.failurestress=[NaN,NaN,stress];
[lpostnotbreak,~,~,xnotbreak]=g_calcprior(ResultSet,[],lprior);
xnotbreak=1-xnotbreak; %x2surv is the cdf of failing at runoutstress, so 
%need to flip to make it survive

%calculate the marginal probability
probbreak=exp(lprior).*xbreak;
probnotbreak=exp(lprior).*xnotbreak;
%integrate
probbreak=sum(probbreak,'all');
probnotbreak=sum(probnotbreak,'all');

%normalise the posterior, and finally the integral across the prior:
normbreak=sum(exp(lpostbreak),'all');
lpostbreak=lpostbreak-ones(size(lprior)).*log(normbreak);
normnotbreak=sum(exp(lpostnotbreak),'all');
lpostnotbreak=lpostnotbreak-ones(size(lprior)).*log(normnotbreak);
%information is p*logp
infobreak=sum(lpostbreak.*exp(lpostbreak),'all');
infonotbreak=sum(lpostnotbreak.*exp(lpostnotbreak),'all');
%output: negative for original convenience of finding the minimum of y
%(which is the maximum of utility)
y=-(probbreak*infobreak+probnotbreak*infonotbreak);

end

