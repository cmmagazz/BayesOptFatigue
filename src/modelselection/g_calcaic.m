function [aic,bic,normloglike,wblloglike,lognloglike]=g_calcaic(data,constantlife1,theta,sigma)
%find loglikelihodd surfaces, AIC and BIC for normal, lognormal, and
%weibull distibutions given the data.

%INPUTS: theta, sigma, data (failurestress or ResultSet), constantlifel (1
%= constant life, 0 = PSN)

%OUTPUTS: AIC, BIC, log likelihood surfaces for normal, weibull,log-normal
%distributions

%wblconsta and lognconst are constants to make prior spaces similar sizes
%compared to parameter distributions

[mtheta,msigma]=meshgrid(theta,sigma);

mtheta=mtheta'; %because meshgrid
msigma=msigma';

wblconst=4;
lognconst=150;

if constantlife1==1
    [normloglike,wblloglike,lognloglike,maxnormloglike,maxwblloglike,maxlognloglike]=g_calcloglike(data,theta,sigma,wblconst,lognconst);
    %DEBUG - plot the log-likelihood surface
    %{
    figure(1)
    subplot(1,3,1)
    h=pcolor(msigma,mtheta,normloglike);
    set(h,'EdgeColor','none')
    caxis([max(normloglike(:))-20,max(normloglike(:))])
    hcb=colorbar;
    ylabel(hcb, 'ln(L)')
%     hcb.Title.String='ln(L)';
    hold on
    [~,id]=max(normloglike(:));
    scatter(msigma(id),mtheta(id),'kx')
%     axis image
    xlabel('\sigma')
    ylabel('\mu')
    title('Normal Distribution')
    subplot(1,3,2)
    h=pcolor(msigma./wblconst,mtheta,wblloglike);
    caxis([max(wblloglike(:))-20,max(wblloglike(:))])
    set(h,'EdgeColor','none')
    hcb=colorbar;
    ylabel(hcb, 'ln(L)')
%     axis image
    hold on
    [~,id]=max(wblloglike(:));
    scatter((msigma(id)./wblconst),mtheta(id),'kx')
    
    
    xlabel('\beta')
    ylabel('\alpha')
    title('Weibull Distribution')
    subplot(1,3,3)
    h=pcolor(msigma./lognconst,log(mtheta),lognloglike);
    hold on
    
    set(h,'EdgeColor','none')
    caxis([max(lognloglike(:))-20,max(lognloglike(:))])
    hcb=colorbar;
    ylabel(hcb, 'ln(L)')
%     axis image
    [~,id]=max(lognloglike(:));
    scatter((msigma(id)./lognconst),log(mtheta(id)),'kx')
    xlabel('\sigma')
    ylabel('\mu')
    title('Log-Normal Distribution')
    drawnow
    %}
    logL=[maxnormloglike;maxwblloglike;maxlognloglike];
    numParam=2;
else
    init=gPSN_calcprior([],theta,sigma);%theta=A, sigma=B
    idx=exp(lprior)==max(exp(lprior));
    logL=lprior(idx)-init(idx);
    numparam=size(A,1)+size(B,1);
end
%find number of observations
try
    numObs=1:size(data.failuretally,1);
catch
    numObs=1:size(data,1);
end
%Calculate AIC
for i=1:size(logL,1)%for every model
    [aic(i,:),bic(i,:)]=aicbic(logL(i,:),numParam,numObs);
end