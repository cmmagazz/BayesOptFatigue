
ResultSet.details=f_setupresultsdist([[150, 650];[1, 400]],'norm',[250,200]);

%% convert old totalresults SN data into new version


totalresults=f_SNresults(totalresults);

%% setup results
bic=NaN(4,6,30);

for i=1:4
    
    if i==1
        [~,bicT]=g_calcaic(staircaseresults.SN.results,1);
    elseif i==2
        [~,bicT]=g_calcaic(bayesresults.SN.results,1);
    elseif i==3
        [~,bicT]=g_calcaic(stepresults.SN.results,1);
    elseif i==4
        [~,bicT]=g_calcaic(bayesstepresults.SN.results,1);
    end
        
        
    bic(i,:,1:size(bicT,2))=bicT;

end
%% plot relative to normal

co='grbcmyrgbcmyk';
figure
for i= 1:4
        subplot(2,2,i)
        bictemp=squeeze(bic(i,:,:));
        k=2;%number of model params
        n=repmat(permute(1:30,[3,1,2]),[1,3,1]);%number of samples for each model

        bictemp=bictemp-repmat(bictemp(3,:),[6,1]);
        plot(bictemp')
        xlabel('Number of samples Tested')
        ylabel('\DeltaBIC')
        if i ==1
            title('Staircase')
        elseif i ==2
            title('Bayesian Staircase')
        elseif i ==3
            title('Stress Step')
        elseif i ==4
            title('Bayesian Stress Step')    
        end
        hold on
        yline(2,'k--');
        yline(6,'k--');
        yline(10,'k--');
        yline(-2,'k-.');
        yline(-6,'k-.');
        yline(-10,'k-.');
        ylim([-12,12])
        sgtitle('Experimental \Delta BIC')

    
end
legend('Normal','Log-Normal','2-Parameter Weibull','3-Parameter Weibull','Gumbell','Frechet')
%%
print(...
    ['experimental bic for each protocol-minusoutlier'], '-dpng','-r1500')
%% for all of it
totalresults.SN.results(95,:)=[];
% totalresults=f_SNresults(totalresults);
        [~,bigbic]=g_calcaic(totalresults.SN.results,1);
        bigbic=bigbic-repmat(bigbic(3,:),[6,1]);

plot(bigbic',LineWidth=1.5)

        xlabel('Number of samples Tested')
        ylabel('\DeltaBIC')

        hold on
        yline(2,'k--');
        yline(6,'k--');
        yline(10,'k--');
        yline(-2,'k-.');
        yline(-6,'k-.');
        yline(-10,'k-.');
        ylim([-12,12])
legend('Normal','Log-Normal','2-Parameter Weibull','Location','SouthWest')%,'3-Parameter Weibull','Gumbell','Frechet')

        %title('Experimental \Delta BIC - All Results')
        
%% save        
print(...
    ['experimental bic for all data - minus sample 95'], '-dpng','-r1500')
savefig('experimental bic for all data - minus sample 95')

