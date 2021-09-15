
ResultSet.details=f_setupresultsdist([[150, 650];[1, 400]],'norm',[250,200]);

%% convert old totalresults SN data into new version


totalresults=f_SNresults(totalresults);

%% setup results
bic=NaN(4,3,30);

for i=1:4
    
    if i==1
        [~,bicT]=g_calcaic(staircaseresults.SN.results,1,ResultSet.details.theta{1},ResultSet.details.theta{2});
    elseif i==2
        [~,bicT]=g_calcaic(bayesresults.SN.results,1,ResultSet.details.theta{1},ResultSet.details.theta{2});
    elseif i==3
        [~,bicT]=g_calcaic(stepresults.SN.results,1,ResultSet.details.theta{1},ResultSet.details.theta{2});
    elseif i==4
        [~,bicT]=g_calcaic(bayesstepresults.SN.results,1,ResultSet.details.theta{1},ResultSet.details.theta{2});
    end
        
        
    bic(i,:,:)=bicT;

end
%% plot relative to normal

co='grbcmyrgbcmyk';
figure
for i= 1:4
        subplot(2,2,i)
        bictemp=squeeze(bic(i,:,:));
        k=2;%number of model params
        n=repmat(permute(1:30,[3,1,2]),[1,3,1]);%number of samples for each model
        correctionterm=(2*k^2+2*k)./(n-k-1);
        bictemp=bictemp-repmat(bictemp(2,:),[3,1]);
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
legend('Normal','Weibull','Log-Normal')
print(...
    ['experimental bic for each protocol'], '-dpng','-r1500')
%% for all of it

% totalresults=f_SNresults(totalresults);
        [~,bigbic]=g_calcaic(totalresults.SN.results,1,ResultSet.details.theta{1},ResultSet.details.theta{2});
        bigbic=bigbic-repmat(bigbic(2,:),[3,1]);
plot(bigbic')
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
        legend('Normal','Weibull','Log-Normal')

        title('Experimental \Delta BIC - All Results')
        
        
print(...
    ['experimental bic for all dta - minus a few'], '-dpng','-r1500')