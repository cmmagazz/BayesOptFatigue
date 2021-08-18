function [loglike,maxloglike]=g_calcloglike(ResultSetraw,theta,distnames,NUMEL)
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
loglike=cell(length(distnames),1);
maxloglike=loglike;
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
        for distid=1:length(distnames)
            dist=distnames{distid};
            ResultSet.details=f_setupresultsdist(theta{distid},dist,NUMEL(distid));
            if ~isnan(failstress(i))&&~isnan(runoutstress(i))%usual case when we have a failure and runout
                x2fail=g_calcprobCDF(failstress(i),ResultSet.details.theta,dist);
                x2surv=g_calcprobCDF(runoutstress(i),ResultSet.details.theta,dist);
                x2=x2fail-x2surv;%the probability of getting sample data in the range is the integral of the
            elseif isnan(runoutstress(i)) %if we failed on the first step
                x2fail=g_calcprobCDF(failstress(i),ResultSet.details.theta,dist);
                x2=x2fail;
            elseif isnan(failstress(i)) %if we reached the limit of testing stresses and did not fail [should not ever happen!]
                x2surv=g_calcprobCDF(runoutstress(i),ResultSet.details.theta,dist);
                x2=1-x2surv;
            end
            %PDF of collection probability over the range. This is equivelant to the difference between 
            %the values of the cdfover that range
            whereuseful=x2>0;
            if i==1
                loglike{distid}=NaN(size(x2));
                loglike{distid}(whereuseful)=log(x2(whereuseful));
            else
                loglike{distid}(whereuseful)=log(x2(whereuseful))+loglike{distid}(whereuseful);
            end
            loglike{distid}(~whereuseful)=-1e7;
            maxloglike{distid}(i)=max(loglike{distid}(:));
        end
    end
else
    error('no data')
end

end