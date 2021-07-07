function [normloglike,wblloglike,lognloglike,maxnormloglike,maxwblloglike,maxlognloglike]=g_calcloglike(ResultSetraw,theta,sigma,wblconst,lognconst)
%Calculate the loglikelihood surfaces for normal, log normal, and weibull
%distributions given the data

%INPUTS: ResultSet.raw (or failurestress),theta,sigma,wblconst 
%(scaling constant),lognconst (scaling constant)
%OUTPUTS: loglikelihood surfacses and maximum values for normal, weibull
%and log-normal distributions

%remove nans and flip failuretally to the engineering definition (pass/fail) definition
if isstruct(ResultSetraw)
    data=ResultSetraw.failurestress(:,:);
else
    data=ResultSetraw;
end
if ~isempty(data)
    idS=isnan(data(:,1)) & isnan(data(:,3));
    data(idS,:)=[];
    %data=f_stripexcessrunouts(data);
end

if ~isempty(data)
    failstress=data(:,1);
    runoutstress=data(:,3);
    %Now create a prior which, if the sample passed, then you have informaiton
    %on the passing cdf (decreases with increased stress), and vice versa for
    %failure. If the probability of this is 1 or 0, set the prior to -100;
    %run through all the data
    %Now create a prior which, if the sample passed, then you have informaiton
    %on the passing cdf (decreases with increased stress), and vice versa for
    %failure. If the probability of this is 1 or 0, set the prior to -100;
    %run through all the data
    for i=1:numel(failstress)
        if ~isnan(failstress(i))&&~isnan(runoutstress(i))%usual case when we have a failure and runout
            normx2fail=normcdf(failstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
            normx2surv=normcdf(runoutstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
            normx2=normx2fail-normx2surv;%the probability of getting sample data in the range is the integral of the
            
            wblx2fail=wblcdf(failstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma./wblconst,[size(theta,2),1]));
            wblx2surv=wblcdf(runoutstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma./wblconst,[size(theta,2),1]));
            wblx2=wblx2fail-wblx2surv;%the probability of getting sample data in the range is the integral of the
            
            lognx2fail=logncdf(failstress(i),repmat(log(theta)',[1,size(sigma,2)]),repmat(sigma./lognconst,[size(theta,2),1]));
            lognx2surv=logncdf(runoutstress(i),repmat(log(theta)',[1,size(sigma,2)]),repmat(sigma./lognconst,[size(theta,2),1]));
            lognx2=lognx2fail-lognx2surv;%the probability of getting sample data in the range is the integral of the
        elseif isnan(runoutstress(i)) %if we failed on the first step
            normx2fail=normcdf(failstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
            normx2=normx2fail;
            
            wblx2fail=wblcdf(failstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma./wblconst,[size(theta,2),1]));
            wblx2=wblx2fail;
            
            lognx2fail=logncdf(failstress(i),repmat(log(theta)',[1,size(sigma,2)]),repmat(sigma./lognconst,[size(theta,2),1]));
            lognx2=lognx2fail;
        elseif isnan(failstress(i)) %if we reached the limit of testing stresses and did not fail [should not ever happen!]
            normx2surv=normcdf(runoutstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
            normx2=1-normx2surv;
            
            wblx2surv=wblcdf(runoutstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma./wblconst,[size(theta,2),1]));
            wblx2=1-wblx2surv;
            
            lognx2surv=logncdf(runoutstress(i),repmat(log(theta)',[1,size(sigma,2)]),repmat(sigma./lognconst,[size(theta,2),1]));
            lognx2=1-lognx2surv;
        end
        %PDF of collection probability over the range. This is equivelant to the difference between 
        %the values of the cdfover that range
        normwhereuseful=normx2>0;
        
        wblwhereuseful=wblx2>0;
        
        lognwhereuseful=lognx2>0;
        if i==1
            normloglike=NaN(size(normx2));
            normloglike(normwhereuseful)=log(normx2(normwhereuseful));
            
            wblloglike=NaN(size(normx2));
            wblloglike(wblwhereuseful)=log(wblx2(wblwhereuseful));
            
            lognloglike=NaN(size(normx2));
            lognloglike(lognwhereuseful)=log(lognx2(lognwhereuseful));
        else
            normloglike(normwhereuseful)=log(normx2(normwhereuseful))+normloglike(normwhereuseful);
            
            wblloglike(wblwhereuseful)=log(wblx2(wblwhereuseful))+wblloglike(wblwhereuseful);
            
            lognloglike(lognwhereuseful)=log(lognx2(lognwhereuseful))+lognloglike(lognwhereuseful);
        end
        normloglike(~normwhereuseful)=-1e7;
        maxnormloglike(i)=max(normloglike(:));
        
        wblloglike(~wblwhereuseful)=-1e7;
        maxwblloglike(i)=max(wblloglike(:));
        
        lognloglike(~lognwhereuseful)=-1e7;
        maxlognloglike(i)=max(lognloglike(:));
    end
else
    error('no data')
end

end