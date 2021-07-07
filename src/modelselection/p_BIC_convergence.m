%Runs many tests on multiple protocols to compare effects on the evidence
%for distributions of the input sample set.
%=============================================
%protocol constants
tic
%protocolstotest={1};
filename='Paper figure BIC with different protocols';
protocolstotest={2,12,5,9};
%protocolstotest={2};
%filename = 'testAnimated3.gif';
ResultSet.details.protocolnames={'adaptive','stress step','random',...
    'adaptive midstep','staircase','probit','life','simple','bayes staircase',...
    'bayes staircase flat prior','bayes staircase 1/x prior','bayes step'};
%setup the results parameters
%=============================================
%step size constants
ResultSet.details.startingstress=280;
ResultSet.details.runout=6; %set the runout value
ResultSet.details.step.steptype=1;%1=abs, 0=percent
ResultSet.details.step.Nstepdown=3;
redefineq=1;
TestingSet.details.SDistq=1; %sample distribution: 1=gaussian, 2=weibull
ResultSet.details.step.stepsize=5;
%=============================================
%SN data constants
widths=[10];%,60];

%=============================================
%prior settings
maxtheta=650;
mintheta=150;
thetaR=mintheta:1:maxtheta;
ResultSet.details.theta=thetaR;
% set upper and lower bound of possible values of sigma here
maxsigma=200;
minsigma=1;
sigmaR=minsigma:1:maxsigma;
ResultSet.details.sigma=sigmaR;
%=================================================
nrepeats=100; %number of times to test
nsamptotest=100;
%nsamptotest=unique(nsamptotest)+1; %equally spaced in log sample numbers, roughly 30 
aic=NaN(length(protocolstotest),numel(widths),nrepeats,3,nsamptotest);
bic=NaN(length(protocolstotest),numel(widths),nrepeats,3,nsamptotest);
sdistk=2;
for protocolk= 1:length(protocolstotest)
    ResultSet.details.protocol = ResultSet.details.protocolnames{protocolstotest{protocolk}};
    if strcmp(ResultSet.details.protocol,'bayes staircase')
        ResultSet.details.startingstress=350;
    elseif strcmp(ResultSet.details.protocol,'staircase')
        ResultSet.details.startingstress=350;
        ResultSet.details.step.stepsize=20;
        
    else
        ResultSet.details.startingstress=100;
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
            TestingSet.details.N=linspace(0.01,10,1000);
            TestingSet.details.DAMq=1;
            TestingSet.details.width=widths(widthl);%./stdconst(sdistk);
            TestingSet.details.numsamp=nsamptotest;
            sigmax=NaN(length(nrepeats),1);
            thetamax=NaN(length(nrepeats),1);
            disp(strcat('protocol ',{' '}, num2str(protocolk),'/',num2str(length(protocolstotest)),' width '...
                ,{' '}, num2str(widthl),'/',num2str(length(widths)),'sample',{' '},' repeat',{' '},...
                num2str(repeatz),'/',num2str(nrepeats)))
            %run the variable step sizes
            TestingSet.meanFS=f_createsample(TestingSet.details);
            %run the protocol
            ResultSetraw=f_testingProtocol(TestingSet,ResultSet);
            [aicT,bicT]=g_calcaic(ResultSetraw,1,thetaR,sigmaR);
            aic(protocolk,widthl,repeatz,:,:)=aicT;
            bic(protocolk,widthl,repeatz,:,:)=bicT;
        end
    end
end
toc

%% plot BIC
co='grbcmyrgbcmyk';
figure
for protocolk= 1:length(protocolstotest)
    for widthl=1%1:numel(widths)
        subplot(2,2,protocolk)
        bictemp=squeeze(bic(protocolk,widthl,:,:,:));
        k=2;%number of model params
        n=repmat(permute(1:nsamptotest,[3,1,2]),[nrepeats,3,1]);%number of samples for each model
        correctionterm=(2*k^2+2*k)./(n-k-1);
        bictemp=bictemp-repmat(bictemp(:,1,:),[1,size(bictemp,2),1]);
        meanarea=squeeze(mean(bictemp,1));
        stdarea=squeeze(std(bictemp,[],1));
        shadedErrorBar(1:nsamptotest,meanarea(1,:),stdarea(1,:),'lineProps',strcat('-',co(1)))
        hold on 
        shadedErrorBar(1:nsamptotest,meanarea(2,:),stdarea(2,:),'lineProps',strcat('-',co(2)))
        shadedErrorBar(1:nsamptotest,meanarea(3,:),stdarea(3,:),'lineProps',strcat('-',co(3)))
        xlabel('Number of samples Tested')
        ylabel('\DeltaBIC')
        if protocolk ==1
            title('Stress Step')
        elseif protocolk ==2
            title('Bayes Step')
        elseif protocolk ==3
            title('Staircase')
        elseif protocolk ==4
            title('Bayes Staircase')    
        end
        hold on
        yline(2,'k--');
        yline(6,'k--');
        yline(10,'k--');
        yline(-2,'k-.');
        yline(-6,'k-.');
        yline(-10,'k-.');
        ylim([-12,12])
        sgtitle(strcat('Weibull Distribution, \beta = ',num2str(widths(widthl))))

    end
end
J=protocolstotest{:}
legend('Normal','Weibull','Log-Normal')

% sgtitle(strcat('Stress Step, Bayes Information Criterion'));%,', Stepsize = ',{' '},num2str(ResultSet.details.step.stepsize)))
% set(gcf, 'Position', get(0, 'Screensize'));
% print([strcat('Number_of_Repeats_',num2str(nrepeats),'normalised AIC convergence, stepsize = ',num2str(ResultSet.details.step.stepsize)) filename], '-dpng','-r300')
% 
% %close all
% 
% currdate=datestr(datetime);
% currdate=currdate(1:11);
% save([strcat('Number_of_Repeats_',num2str(nrepeats),'convergence test  50 test oone width') filename currdate '.mat']);
