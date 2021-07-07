function [y,secondoutput]= g_STEP_UTILITY(stepsize,startingstress,theta,sigma,lprior)
% Utility function based on Gearhart method
% Conceptually: we want to find a stress and stress step where, whether it breaks or fails (at either stress),
% it will give us the highest information gain. This is also scaled by how
% likely it is, given our previous data, for a pass or fail. 
%
% input:
%   stepsize: a double 
%   startingstress: a double
%   theta: space of theta values allowed, vector
%   sigma: space of sigma values allowed, vector
%   lprior: prior distribution across the theta/sigma space
%
% output:
%   y: value for utility at this stress, given this prior distribution
%
% CMM RR 2021
priorinfo=sum(lprior.*exp(lprior),'all');
minstress=startingstress;
maxval=f_HPD(squeeze(lprior),1-1e-6);
[row,col]=find(squeeze(exp(lprior))>=maxval);
thetarange=theta(row);
sigmarange=sigma(col);
maxstress=max(thetarange)+max(sigmarange)+2*stepsize+100;
minstressT=min(thetarange)-max(sigmarange)-2*stepsize-100;
%deal with max and min stresses outside a 'reasonable' range
if isnan(minstressT)
    warning('error in calculation')
end
if minstressT<minstress
    1;
else
%     round the min stress values onto the 'grid'
    minstress=interp1(minstress-100*stepsize:stepsize:max(theta),minstress-100*stepsize:stepsize:max(theta),minstressT,'nearest');
    %we only want to look 
%     within the 95% HPD of the prior, but since we want to 
%     controll our pertubation of starting stress about the 'grid' of testing 
%     values to be able to extract more information, we want to round the minimum 
%     of the prior towards values that we would have tested on the perturbed 'grid'
end

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
        [~,failnormstcdf,infobreakT,~]=g_calcprior([stresses(ids),NaN,NaN],theta,sigma,lprior);
        probtotT=squeeze(exp(lprior)).*failnormstcdf;
        probtot(1)=sum(probtotT(:));
    else
%         probability of oberving is the probability of recording a result
%         between the stress level and the stress level + step size
        [~,failnormstcdf,infobreakT,runoutnormstcdf]=g_calcprior([stresses(ids),NaN,stresses(ids)-stepsize],theta,sigma,lprior);
        xtot=failnormstcdf-runoutnormstcdf;
        probtotT=squeeze(exp(lprior)).*xtot;
        probtot(ids)=sum(probtotT(:));
    end
    infobreak(ids)=infobreakT-priorinfo;%information is now defined as information geined by the test to make sum work properly
    infobreak(probtot==0)=0;%deal with NaN*0=inf problems
    
end
UtilityS=(infobreak.*probtot);
%output:
y=-1*sum(UtilityS);
secondoutput=1;%infobreak(3);
end
