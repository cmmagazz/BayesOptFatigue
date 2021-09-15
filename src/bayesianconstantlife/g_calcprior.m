function [lprior,x2fail,shannon,x2surv,progthet,progsig]=g_calcprior(ResultSet,theta,varargin)
% Calculate the prior, or joint posterior.
% Here we initialise, and evaluate the contribution of each sample onto the
% joint posterior given the results. Assumes constant step size, should
% there be step testing.
%
% INPUTS
%   ResultSetraw: struct with only the raw field of ResultSet, OR the
%   failurestress table immediately. 
%   theta: the values possible, as a vector
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
if isstruct(ResultSet)
    data=ResultSet.raw.failurestress(:,:);
    theta=ResultSet.details.theta;
    dist=ResultSet.details.dist;
else
    data=ResultSet;
    dist='norm';
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
        psigma=1./(theta{2});
        psigma=psigma./sum(psigma);
        pmu=ones(size(theta{1}))./length(theta{1});
        lprior=log(psigma.*pmu');
    elseif lpriorq==2
        if numel(theta)==2
            lprior=ones(length(theta{1}),length(theta{2}));
        elseif numel(theta)==3
            lprior=ones(length(theta{1}),length(theta{2}),length(theta{3}));
        end
    end

    norm=sum(sum(exp(lprior)));
    lprior=lprior-ones(size(lprior)).*log(norm);
    shannon=sum(lprior.*exp(lprior),'all');
    x2fail=[];
    x2surv=[];
else
    if ~isempty(data)
        data=data(end,:);
    end
end
    
if ~isempty(data)
    failstress=data(:,1);
    runoutstress=data(:,3);
    %Now create a prior which, if the sample passed, then you have informaiton
    %on the passing cdf (decreases with increased stress), and vice versa for
    %failure. If the probability of this is 1 or 0, set the prior to a
    %very small number
    %run through all the data

    shannon = NaN(1,numel(failstress));

    progthet=NaN(numel(failstress),3);
    progsig=NaN(numel(failstress),3);
    progC=NaN(numel(failstress),3);
    for i=1:numel(failstress)
        if ~isnan(failstress(i))&&~isnan(runoutstress(i))%usual case when we have a failure and runout
            x2fail=g_calcprobCDF(failstress(i),theta,dist);
            x2surv=g_calcprobCDF(runoutstress(i),theta,dist);
            x2=x2fail-x2surv;%the probability of getting sample data in the range is the integral of the 
        elseif isnan(runoutstress(i)) %if we failed on the first step
            x2fail=g_calcprobCDF(failstress(i),theta,dist);
            x2=x2fail;
            x2surv=[];
        elseif isnan(failstress(i)) %if we reached the limit of testing stresses and did not fail 
            x2surv=g_calcprobCDF(runoutstress(i),theta,dist);
            x2=1-x2surv;
            x2fail=[];
        end
        %PDF of collection probability over the range. This is equivelant to the difference between 
        %the values of the cdfover that range
        whereuseful=x2>0;
        lprior(whereuseful)=log(x2(whereuseful))+lprior(whereuseful);
        lprior(~whereuseful)=-1e7;

        norm=sum(exp(lprior),'all');
        lprior=lprior-ones(size(lprior)).*log(norm);
        shannon(i)=sum(lprior.*exp(lprior),'all');

        %DEBUG
        %{
        figure(2)
        subplot(2,1,1)
        h=pcolor(theta{2},theta{1},x2);
        set(h, 'EdgeColor', 'none');
        colormap(parula)
        xlabel('\sigma')
        ylabel('\theta')
        title('X2')
        subplot(2,1,2)
        h=pcolor(theta{2},theta{1},exp(lprior));
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
                %find the maximum extent, in theta and sigma independantly,
        %that this draws
        %save the data:
        if ndims(idmin)==2
            [row,col]=find(idmin);

            progthet(i,2)=max(theta{1}(row));
            progthet(i,3)=min(theta{1}(row));
            progsig(i,2)=max(theta{2}(col));
            progsig(i,3)=min(theta{2}(col));
                    %find the theta and sigma with the highest value in lprior
            maxval=max(exp(lprior(:))); %<<<<<<<<< CMM EDIT FROM 17/06/21 this is right!
            idmin=exp(lprior)==maxval;
            [row,col]=find(idmin);
            if nnz(row)>0
                progthet(i,1)=mean(theta{1}(row));
                progsig(i,1)=mean(theta{2}(col));
            end
        elseif ndims(idmin)==3
            I=find(idmin);
            [row,col,depth]=ind2sub(size(idmin),I);
            progthet(i,2)=max(theta{1}(row));
            progthet(i,3)=min(theta{1}(row));
            progsig(i,2)=max(theta{2}(col));
            progsig(i,3)=min(theta{2}(col));
            progC(i,2)=max(theta{3}(depth));
            progC(i,3)=min(theta{3}(depth));
                    %find the theta and sigma with the highest value in lprior
            maxval=max(exp(lprior(:))); %<<<<<<<<< CMM EDIT FROM 17/06/21 this is right!
            idmin=exp(lprior)==maxval;
            
            I=find(idmin);
            [row,col,depth]=ind2sub(size(idmin),I);

            if nnz(row)>0
                progthet(i,1)=mean(theta{1}(row));
                progsig(i,1)=mean(theta{2}(col));
                progC(i,3)=mean(theta{3}(depth));
            end
        end



    end
    norm=sum(exp(lprior),'all');
    lprior=lprior-ones(size(lprior)).*log(norm);
else
    x2fail=[];
    x2surv=[];
    progthet=[];
    progsig=[];
end
    
norm=sum(exp(lprior),'all'); %Sanity: normalise at the end as well
lprior=lprior-ones(size(lprior)).*log(norm);
if nnz(isnan(lprior))>0 
    1;
    warning('Nans in the prior!')
end
end