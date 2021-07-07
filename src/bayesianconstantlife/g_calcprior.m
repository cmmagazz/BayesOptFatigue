function [lprior,x2fail,shannon,x2surv,progthet,progsig]=g_calcprior(ResultSetraw,theta,sigma,varargin)
% Calculate the prior, or joint posterior.
% Here we initialise, and evaluate the contribution of each sample onto the
% joint posterior given the results. Assumes constant step size, should
% there be step testing.
%
% INPUTS
%   ResultSetraw: struct with only the raw field of ResultSet, OR the
%   failurestress table immediately. 
%   theta: the mean values possible, as a vector
%   sigma: the standard deviation values possible, as a vector
%   varargin: 
%       lprior: the log of the prior/joint posterior, as an array. if not
%       given, we initialise and evaluate
%       lpriorq: a flag for how to initialise:
%               1 = Jeffreys or 1/sigma
%               2 = flat
%
% OUTPUTS
%       lprior: the log of the prior, as an array
%       x2fail: the likelihood of recording a failure given the parameters
%       shannon: the shannon information (=p*log(p))
%       x2surv: 1-the likelihood of recording a runout given the parameters
%       progthet: the progressive estimation of the MLE,95%CI of the mean
%       progsig: the progressive estimation of the MLE,95%CI of the stdev

%Parse the inputs:
if isstruct(ResultSetraw)
    data=ResultSetraw.failurestress(:,:);
else
    data=ResultSetraw;
end
if ~isempty(data)
    idS=isnan(data(:,1)) & isnan(data(:,3));
    data(idS,:)=[];
end
if numel(varargin)==2
    if size(varargin{1},1)==1
        lpriorq=varargin{1};
        lprior=varargin{2};
    else
        lprior=varargin{1};
        lpriorq=varargin{2};
    end
elseif numel(varargin)==1
    if size(varargin{1},1)==1
        lpriorq=varargin{1};
    else
        lprior=varargin{1};
    end
end
if ~exist('lprior','var') %if you are not given a prior:do the old method where you reupdate
    if ~exist('lpriorq','var')%set deafult to what we have always done
        lpriorq=2;
    end
    if lpriorq==1
        psigma=1./(sigma);
        psigma=psigma./sum(psigma);
        ptheta=ones(size(theta))./length(theta);
        lprior=log(psigma.*ptheta');
    elseif lpriorq==2
        lprior=ones(length(theta),length(sigma));
    end

    norm=sum(sum(exp(lprior)));
    lprior=lprior-ones(size(lprior)).*log(norm);
    shannon=sum(lprior.*exp(lprior),'all');
    x2fail=[];
    x2surv=[];
    if ~isempty(data)
        failstress=data(:,1);
        runoutstress=data(:,3);
        %Now create a prior which, if the sample passed, then you have informaiton
        %on the passing cdf (decreases with increased stress), and vice versa for
        %failure. If the probability of this is 1 or 0, set the prior to a
        %very small number
        %run through all the data
        shannonT=shannon;
        shannon = NaN(1,numel(failstress));
        shannon(1)=shannonT;
        progthet=NaN(numel(failstress),3);
        progsig=NaN(numel(failstress),3);
        for i=1:numel(failstress)
            if ~isnan(failstress(i))&&~isnan(runoutstress(i))%usual case when we have a failure and runout
                x2fail=normcdf(failstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
                x2surv=normcdf(runoutstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
                x2=x2fail-x2surv;%the probability of getting sample data in the range is the integral of the 
            elseif isnan(runoutstress(i)) %if we failed on the first step
                x2fail=normcdf(failstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
                x2=x2fail;
            elseif isnan(failstress(i)) %if we reached the limit of testing stresses and did not fail 
                x2surv=normcdf(runoutstress(i),repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
                x2=1-x2surv;
            end
            %PDF of collection probability over the range. This is equivelant to the difference between 
            %the values of the cdfover that range
            whereuseful=x2>0;
            lprior(whereuseful)=log(x2(whereuseful))+lprior(whereuseful);
            lprior(~whereuseful)=-1e7;
            
            norm=sum(exp(lprior(:)));
            lprior=lprior-ones(size(lprior)).*log(norm);
            shannon(i+1)=sum(lprior.*exp(lprior),'all');
            
            %DEBUG
            %{
            figure(2)
            subplot(2,1,1)
            h=pcolor(sigma,theta,x2);
            set(h, 'EdgeColor', 'none');
            colormap(parula)
            xlabel('\sigma')
            ylabel('\theta')
            title('X2')
            subplot(2,1,2)
            h=pcolor(sigma,theta,exp(lprior));
            set(h, 'EdgeColor', 'none');
            colormap(parula)
            xlabel('\sigma')
            ylabel('\theta')
            xlim([0,80])
            ylim([300,500])
            title('lprior')
            sgtitle('g calcprior')
            pause(0.5)
            %}
            
            %Given the full dataset, can also output the progressive
            %estimation of the MLE and 95% confidence interval as it
            %evaluates
            CI=0.95;
            maxval=f_HPD(lprior,CI);
            idmin=exp(lprior)>=maxval;
            [row,col]=find(idmin);
            %find the maximum extent, in theta and sigma independantly,
            %that this draws
            %save the data:
            progthet(i,2)=max(theta(row));
            progthet(i,3)=min(theta(row));
            progsig(i,2)=max(sigma(col));
            progsig(i,3)=min(sigma(col));
            %find the theta and sigma with the highest value in lprior
            maxval=max(exp(lprior(:))); %<<<<<<<<< CMM EDIT FROM 17/06/21 this is right!
            idmin=exp(lprior)==maxval;
            [row,col]=find(idmin);
            if nnz(row)>0
                progthet(i,1)=mean(theta(row));
                progsig(i,1)=mean(sigma(col));
            end
            
        end
        norm=sum(exp(lprior(:)));
        lprior=lprior-ones(size(lprior)).*log(norm);
    end
    
else %otherwise, if given a prior, then it must be already with the data
    %from the previous tests. so, just take the last one
    shannon(1)=sum(lprior.*exp(lprior),'all');
    x2fail=[];
    x2surv=[];
    if ~isempty(data)
        failstress=data(end,1);
        runoutstress=data(end,3);
        if ~isnan(failstress)&&~isnan(runoutstress)%usual case when we have a failure and runout
            x2fail=normcdf(failstress,repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
            x2surv=normcdf(runoutstress,repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
            x2=x2fail-x2surv;%the probability of getting sample data in the range is the integral of the 
        elseif isnan(runoutstress) %if we failed on the first step
            x2fail=normcdf(failstress,repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
            x2=x2fail;
        elseif isnan(failstress) %if we reached the limit of testing stresses and did not fail 
            x2surv=normcdf(runoutstress,repmat(theta',[1,size(sigma,2)]),repmat(sigma,[size(theta,2),1]));
            x2=1-x2surv;
        end
        %PDF of collection probability over the range. This is equivelant to the difference between 
        %the values of the cdfover that range
        whereuseful=x2>0;
        lprior(whereuseful)=log(x2(whereuseful))+lprior(whereuseful);
        lprior(~whereuseful)=-1e7;
        norm=sum(sum(exp(lprior)));
        lprior=lprior-ones(size(lprior)).*log(norm);
        shannon(2)=sum(lprior.*exp(lprior),'all');
    end
    shannon=shannon(end);
end
norm=sum(sum(exp(lprior))); %Sanity: normalise at the end as well
lprior=lprior-ones(size(lprior)).*log(norm);
if nnz(isnan(lprior))>0 
    1;
    warning('Nans in the prior!')
end
end