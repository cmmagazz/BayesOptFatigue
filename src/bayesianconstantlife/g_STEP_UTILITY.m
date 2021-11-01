function [y,secondoutput]= g_STEP_UTILITY(startingstress,ResultSet)
% Utility function based on Gearhart method
% Conceptually: we want to find a stress and stress step where, whether it breaks or fails (at either stress),
% it will give us the highest information gain. This is also scaled by how
% likely it is, given our previous data, for a pass or fail. 
%
% input:
%   stepsize: a double 
%   startingstress: a double
%   theta: space of parameter values allowed, cell array of vectors
%   lprior: prior distribution across the theta/sigma space
%
% output:
%   y: value for utility at this stress, given this prior distribution
%
% CMM RR 2021
theta=ResultSet.details.theta;
lprior=ResultSet.raw.lprior;
stepsize=ResultSet.details.step.stepsize;

priorinfo=sum(lprior.*exp(lprior),'all');

minstress=startingstress;

murange=theta{1}(:);
maxstress=max(murange)+6*max(theta{2}(:));

stresses=minstress:stepsize:maxstress;

probtot=zeros(length(stresses),1);
infobreak=zeros(length(stresses),1);
for ids=1:length(stresses)
    %We want to know the probability that the sample will fail
    %and if it does fail, what is the information that we get from that failure in
    %the prior
    %if it runs out, we have no extra information as we can't record this
    %with the step testing protocol
    if ids==1%on the first step, the probability of failure is the integral from zero to the starting stress
        ResultSet.raw.failurestress=[stresses(ids),NaN,NaN];
        [~,failnormstcdf,infobreakT,~]=g_calcprior(ResultSet,[],lprior);
        probtotT=squeeze(exp(lprior)).*failnormstcdf;
        probtot(1)=sum(probtotT,'all');
    else
%         probability of oberving is the probability of recording a result
%         between the stress level and the stress level + step size
        ResultSet.raw.failurestress=[stresses(ids),NaN,stresses(ids)-stepsize];

        [~,failnormstcdf,infobreakT,runoutnormstcdf]=g_calcprior(ResultSet,[],lprior);
        xtot=failnormstcdf-runoutnormstcdf;
        probtotT=squeeze(exp(lprior)).*xtot;
        probtot(ids)=sum(probtotT,'all');
    end
    infobreak(ids)=infobreakT-priorinfo;%information is now defined as information geined by the test to make sum work properly
    infobreak(probtot==0)=0;%deal with NaN*0=inf problems
    
end
UtilityS=(infobreak.*probtot);
%output:
y=-1*sum(UtilityS);
secondoutput=1;%infobreak(3);
end
