function [aic,bic,normloglike,wblloglike,lognloglike]=g_calcaic(data,constantlife1)
%find loglikelihodd surfaces, AIC and BIC for normal, lognormal, and
%weibull distibutions given the data.

%INPUTS: theta, sigma, data (failurestress or ResultSet), constantlifel (1
%= constant life, 0 = PSN)

%OUTPUTS: AIC, BIC, log likelihood surfaces for normal, weibull,log-normal
%distributions

%wblconsta and lognconst are constants to make prior spaces similar sizes
%compared to parameter distributions



wblconst=1;
lognconst=1;
distnames={'norm','lognorm','2pwbl','3pwbl','gev','type1'};
theta={[[150, 650];[1, 150]],...
    [[150, 650];[0, 0.5]./lognconst],...
    [[150, 650];[1, 80]./wblconst],...
    [[150, 2000];[1, 150]./wblconst;[-500,500]],...
    [[150, 650];[1, 150]./wblconst;[-5,5]],...
    [[150, 650];[1, 150]]};
NUMEL=[100,100,100,100,100,100];
numparam=[2,2,2,3,3,2];
numplots=numparam-1;
if constantlife1==1
    [loglike,logL]=g_calcloglike(data,theta,distnames,NUMEL);
    %DEBUG - plot the log-likelihood surface
    %{
    figure(1)
    for distid=1:length(distnames)  
        TempResultdetails=f_setupresultsdist(theta{distid},distnames{distid},NUMEL(distid));
        
        if numparam(distid)==2
            subplot(1,sum(numplots),sum(numplots(1:distid)))
            [mtheta,msigma]=meshgrid(TempResultdetails.theta{1},TempResultdetails.theta{2});
            mtheta=mtheta'; %because meshgrid
            msigma=msigma';

            h=pcolor(msigma,mtheta,loglike{distid});
            set(h,'EdgeColor','none')
            caxis([max(loglike{distid}(:))-20,max(loglike{distid}(:))])
            hcb=colorbar;
            ylabel(hcb, 'ln(L)')
        %     hcb.Title.String='ln(L)';
            hold on
            [~,id]=max(loglike{distid}(:));
            scatter(msigma(id),mtheta(id),'kx')
        %     axis image
            xlabel('\sigma')
            ylabel('\mu')
        elseif numparam(distid)==3
            subplot(1,sum(numplots),sum(numplots(1:distid))-1)
            [A,B,C]=meshgrid(TempResultdetails.theta{1},TempResultdetails.theta{2},TempResultdetails.theta{3});
            A=permute(A,[2,1,3]);
            B=permute(B,[2,1,3]);
            C=permute(C,[2,1,3]);
            plotarray1=log(sum(exp(loglike{distid}),3));
            h=pcolor(A(:,:,1),B(:,:,1),plotarray1);
            set(h,'EdgeColor','none')
            caxis([max(loglike{distid}(:))-20,max(loglike{distid}(:))])
            hcb=colorbar;
            ylabel(hcb, 'ln(L)')
        %     hcb.Title.String='ln(L)';
            hold on
            subplot(1,sum(numplots),sum(numplots(1:distid)))
            plotarray2=log(sum(exp(loglike{distid}),2));
            h=pcolor(squeeze(A(:,1,:)),squeeze(C(:,1,:)),squeeze(plotarray2));
            set(h,'EdgeColor','none')
            caxis([max(loglike{distid}(:))-20,max(loglike{distid}(:))])
            hcb=colorbar;
            ylabel(hcb, 'ln(L)')
        %     hcb.Title.String='ln(L)';
            hold on
        end
        title([distnames{distid},' distribution'])
    end
    sgtitle('Log L comparison')
    drawnow
    %}
    %% 
    %{
    figure(2)
    for distid=1:length(distnames)  
        TempResultdetails=f_setupresultsdist(theta{distid},distnames{distid},NUMEL(distid));
        [~,fitid]=max(loglike{distid},[],'all','linear');
        [A,B,C]=ind2sub(size(loglike{distid}),fitid);
        try
            CDF=g_calcprobCDF(0:1000,{TempResultdetails.theta{1}(A),TempResultdetails.theta{2}(B),TempResultdetails.theta{3}(C)},distnames{distid});
        catch
            CDF=g_calcprobCDF(0:1000,{TempResultdetails.theta{1}(A),TempResultdetails.theta{2}(B)},distnames{distid});
        end
        plot(0:1000,gradient(CDF))
        hold on
    end
    legend(distnames)
    drawnow
    ylabel('Probability')
    xlabel('Stress')
    %}
    %{
    figure(3)
    for distid=1:length(distnames)
        subplot(1,length(distnames),distid)
        TempResultdetails=f_setupresultsdist(theta{distid},distnames{distid},NUMEL(distid));
        fitid=find(exp(loglike{distid})>(max(exp(loglike{distid}),[],'all')*0.2));
        [A,B,C]=ind2sub(size(loglike{distid}),fitid);
        for i=1:length(A)
            try
                CDF=g_calcprobCDF(0:1000,{TempResultdetails.theta{1}(A(i)),TempResultdetails.theta{2}(B(i)),TempResultdetails.theta{3}(C(i))},distnames{distid});
            catch
                CDF=g_calcprobCDF(0:1000,{TempResultdetails.theta{1}(A(i)),TempResultdetails.theta{2}(B(i))},distnames{distid});
            end
            plot(0:1000,gradient(CDF))
            hold on
        end
        title(distnames{distid})
%         xlim([300,550])
        ylabel('Probability')
        xlabel('Stress')
    end
    
    %}
else
    init=gPSN_calcprior([],theta,sigma);%theta=A, sigma=B
    idx=exp(lprior)==max(exp(lprior));
    logL=lprior(idx)-init(idx);
    numparam=size(A,1)+size(B,1);
end
%find number of observations
try
    numObs=1:size(data.failurestress,1);
catch
    numObs=1:size(data,1);
end
%Calculate AIC
for distid=1:length(distnames)%for every model
    bic(distid,:)=numparam(distid).*log(numObs)-2.*logL{distid};
    aic(distid,:)=2.*numparam(distid)-2.*logL{distid};

end
%% DEBUG - plot BICs
%{
figure(4)
refdist='norm';
for distid=1:length(distnames)
    bictemp=bic(distid,:)-bic(strcmp(distnames,refdist),:);
    plot(bictemp)
    hold on
end
legend(distnames,'AutoUpdate','off')
xlabel('Number of Samples Tested')
ylabel('\Delta BIC')
yline(2,'k--');
yline(6,'k--');
yline(10,'k--');
yline(-2,'k-.');
yline(-6,'k-.');
yline(-10,'k-.');
ylim([-12,12])
%}