%Runs many tests on multiple protocols to compare effects on the evidence
%for distributions of the input sample set.
%=============================================
%protocol constants
clear
tic
%protocolstotest={1};
filename='Paper figure BIC with different protocols';
protocolstotest={5,9,2,12};
%protocolstotest={2};
%setup the results parameters
%=============================================
%step size constants
redefineq=1;
%=============================================
%SN data constants
widths=[40];%,60];
%=============================================
%prior settings
ResultSet.details=f_setupresultsdist([[150, 650];[1, 400]],'norm',[250,200]);
ResultSet.details.runout=6; %set the runout value

ResultSet.details.protocolnames={'adaptive','stress step','random',...
    'adaptive midstep','staircase','probit','life','simple','bayes staircase',...
    'bayes staircase flat prior','bayes staircase 1/x prior','bayes step'};
%=================================================
nrepeats=100; %number of times to test
nsamptotest=80;
%nsamptotest=unique(nsamptotest)+1; %equally spaced in log sample numbers, roughly 30 
aic=NaN(length(protocolstotest),numel(widths),nrepeats,6,nsamptotest);
bic=NaN(length(protocolstotest),numel(widths),nrepeats,6,nsamptotest);
sdistk=1;
for protocolk=1:length(protocolstotest)
    ResultSet.details.protocol = ResultSet.details.protocolnames{protocolstotest{protocolk}};
    if strcmp(ResultSet.details.protocol,'bayes staircase')
        ResultSet.details.startingstress=380;
    elseif strcmp(ResultSet.details.protocol,'staircase')
        ResultSet.details.startingstress=380;
        ResultSet.details.step.stepsize=40;
    else
        ResultSet.details.startingstress=0;
        ResultSet.details.step.stepsize=40;
    end
    stdconst=[1,4,350];
    for widthl=1:numel(widths)
        parfor repeatz=1:nrepeats
            TestingSet=struct();
            TestingSet.details.SDistq=sdistk; %sample distribution: 1=gaussian, 2=weibull
            %=============================================
            %SN data constants
            TestingSet.details.basquin.c=300;
            TestingSet.details.basquin.alpha=600;
            TestingSet.details.basquin.beta=-log(1/6)/log(6);
            TestingSet.details.N=linspace(1,10,1000);
            TestingSet.details.DAMq=1;
            TestingSet.details.width=widths(widthl);%./stdconst(sdistk);
            TestingSet.details.numsamp=nsamptotest;
%             sigmax=NaN(length(nrepeats),1);
%             thetamax=NaN(length(nrepeats),1);
            disp(strcat('protocol ',{' '}, num2str(protocolk),'/',num2str(length(protocolstotest)),' width '...
                ,{' '}, num2str(widthl),'/',num2str(length(widths)),'sample',{' '},' repeat',{' '},...
                num2str(repeatz),'/',num2str(nrepeats)))
            %run the variable step sizes
            TestingSet.meanFS=f_createsample(TestingSet.details);
            %run the protocol
            ResultSet2=f_testingProtocol(TestingSet,ResultSet);
            [aicT,bicT]=g_calcaic(ResultSet2.raw,1);
            aic(protocolk,widthl,repeatz,:,:)=aicT;
            bic(protocolk,widthl,repeatz,:,:)=bicT;
        end
    end
end
toc

%% plot BIC
colours=['grbcmyrgbcmyk']';


figure
for protocolk= 1:length(protocolstotest)
    for widthl=1%1:numel(widths)
        subplot(2,2,protocolk)
        bictemp=squeeze(bic(protocolk,widthl,:,:,:));
        k=2;%number of model params
        n=repmat(permute(1:nsamptotest,[3,1,2]),[nrepeats,3,1]);%number of samples for each model
        bictemp=bictemp-repmat(bictemp(:,3,:),[1,size(bictemp,2),1]);
        meanarea=squeeze(mean(bictemp,1));
        stdarea=squeeze(std(bictemp,[],1));

        shadedErrorBar(1:nsamptotest,meanarea(1,:),stdarea(1,:),'lineProps',...
            {strcat('-',colours(1,:)),'markerfacecolor',colours(1,:)})
        hold on 
        shadedErrorBar(1:nsamptotest,meanarea(2,:),stdarea(2,:),'lineProps',...
            {strcat('-',colours(2,:)),'markerfacecolor',colours(2,:)})
        shadedErrorBar(1:nsamptotest,meanarea(3,:),stdarea(3,:),'lineProps',...
            {strcat('-',colours(3,:)),'markerfacecolor',colours(3,:)})
%         shadedErrorBar(1:nsamptotest,meanarea(4,:),stdarea(4,:),'lineProps',...
%             {strcat('-',colours(4,:)),'markerfacecolor',colours(4,:)})
%         shadedErrorBar(1:nsamptotest,meanarea(5,:),stdarea(5,:),'lineProps',...
%             {strcat('-',colours(5,:)),'markerfacecolor',colours(5,:)})
%         shadedErrorBar(1:nsamptotest,meanarea(6,:),stdarea(6,:),'lineProps',...
%             {strcat('-',colours(6,:)),'markerfacecolor',colours(6,:)})
        xlabel('Number of samples Tested')
        ylabel('\DeltaBIC')
        if protocolk ==1
            title('Staircase')
        elseif protocolk ==2
            title('Bayesian Staircase')
        elseif protocolk ==3
            title('Stress Step')
        elseif protocolk ==4
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
        sgtitle(strcat('Normal Distribution, \sigma = ',num2str(widths(widthl))))

    end
end

% legend('Normal','Log-Normal','2-Parameter Weibull','3-Parameter Weibull','Gumbell','Frechet')
legend('Normal','Log-Normal','Weibull');%,'3-Parameter Weibull','Gumbell','Frechet')

% sgtitle(strcat('Stress Step, Bayes Information Criterion'));%,', Stepsize = ',{' '},num2str(ResultSet.details.step.stepsize)))
% set(gcf, 'Position', get(0, 'Screensize'));
% print([strcat('Number_of_Repeats_',num2str(nrepeats),'normalised AIC convergence, stepsize = ',num2str(ResultSet.details.step.stepsize)) filename], '-dpng','-r300')
print(...
    [strcat('Number_of_Repeats_',num2str(nrepeats),'normalised BIC convergence, stepsize = '...
    ,num2str(ResultSet.details.step.stepsize)) filename], '-dpng','-r1500')
%%
close all
% 
currdate=datestr(datetime);
currdate=currdate(1:11);
save([strcat('Number_of_Repeats_',num2str(nrepeats),'convergence test  50 test oone width') filename currdate '.mat']);
