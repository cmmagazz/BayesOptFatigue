%% Convergence scripts
% Scripts to reproduce the figures from the paper to compare the error
% after N samples. 
%=============================================
%protocol constants
clear; clc; close all

protocolstotest={1,6};

% 1=stress step
% 2=staircase
% 3=probit
% 4=life
% 5=bayes staircase
% 6=bayes step

%=================================================
nrepeats=200; %number of repeats of the experiment
numsamp=30;
distwidth=40; %true distribution width

%what ratios of width to step size to test
% widths=linspace(0.5,15,10);%Staircase or
widths=[0.5,5,10];%Step


%Set up results arrays
progthet=zeros(length(protocolstotest),numel(widths),nrepeats,numsamp,3);
progsig=zeros(length(protocolstotest),numel(widths),nrepeats,numsamp,3);
shannon=zeros(length(protocolstotest),numel(widths),nrepeats,numsamp);
bigresultstemp=cell(numel(widths),nrepeats);
bigresults=cell(length(protocolstotest));
%run through each protocol
for protocolk=1:length(protocolstotest)
    disp(['Doing protocol: ', num2str(protocolk)])
    progthetT=zeros(numel(widths),nrepeats,numsamp,3);
    progsigT=zeros(numel(widths),nrepeats,numsamp,3);
    %run through each width:step size ratio
    parfor widthl=1:numel(widths)
        %setup testingset
        TestingSet=struct();
        TestingSet.details.SDistq=1; %sample distribution: 1=gaussian, 2=weibull
        TestingSet.details.basquin.c=300;
        TestingSet.details.basquin.alpha=600;
        TestingSet.details.basquin.beta=-log(1/6)/log(6);
        TestingSet.details.N=linspace(0.01,10,100);
        TestingSet.details.DAMq=1;
        TestingSet.details.width=distwidth;
        TestingSet.details.numsamp=numsamp;
        %setup resultset
        ResultSetdetails = struct();
        ResultSetdetails.protocolnames={'stress step',...
            'staircase','probit','life','bayes staircase','bayes step'};
        ResultSetdetails.protocol = ResultSetdetails.protocolnames{protocolstotest{protocolk}};
        ResultSetdetails.runout=6; %set the runout value
        ResultSetdetails.dist='norm';

        maxtheta=1000;
        mintheta=-101;
        theta=mintheta:5:maxtheta;
        % set upper and lower bound of possible values of sigma here
        maxsigma=200;
        minsigma=1;
        sigma=minsigma:4:maxsigma;
        ResultSetdetails.theta{1}=theta;
        ResultSetdetails.theta{2}=sigma;
        
        %CHANGE THE MAIN VARIABLE
        ResultSetdetails.step.stepsize=widths(widthl)*TestingSet.details.width;
        
        %run through each repeat experiment
        for repeatz=1:nrepeats
            disp(['repeat ',num2str(repeatz)])
            %randomise the starting point for step a little
            ResultSetdetails.startingstress=50-rand*ResultSetdetails.step.stepsize;
            %randomise the starting point for staircase a little
            if strcmp(ResultSetdetails.protocol,'staircase') 
                ResultSetdetails.startingstress=405-2*(0.5-rand)*ResultSetdetails.step.stepsize;
            end
            %create a new sample set each time
            TestingSet.meanFS=f_createsample(TestingSet.details);
            %run the protocol
            tempoj1=ResultSetdetails;
            tempoj1.details=tempoj1;
            ResultSettemp=f_testingProtocol(TestingSet,tempoj1);
            ResultSetdetails.raw=ResultSettemp.raw;
            %output some useful things
            [~,~,shannonT(widthl,repeatz,:),~,progthetT(widthl,repeatz,:,:),progsigT(widthl,repeatz,:,:)]=...
                g_calcprior(ResultSetdetails.raw.failurestress,ResultSetdetails.theta,2);
            bigresultstemp{widthl,repeatz}=ResultSetdetails.raw.failurestress;
        end
        disp(['done width ',num2str(widthl)])
    end
    bigresults{protocolk}=bigresultstemp;
    progthet(protocolk,:,:,:,:)=progthetT;
    progsig(protocolk,:,:,:,:)=progsigT;
%     shannon(protocolk,:,:,:)=shannonT(:,:,2:end);
end
disp('Convergence test finished')

save('convergenceredo_200rep_30samp_largeprior.mat')
 %% Plot the figures
 % COMMENT/UNCOMMENT AS NEEDED FOR STEP vs STAIR METHODS
 
 %Define convergence by what's the error after N samples in mean and stdev 

%error in theta vs sigma: relative to the true value which is 400MPa +/-
%distwidth
errbarthet=(squeeze(progthet(:,:,:,:,1))-400).*(100/400);
errbarsig=(squeeze(progsig(:,:,:,:,1))-distwidth).*(100/distwidth);

%START WITH THETA / MEAN
arraytotest=sqrt(errbarthet(:,:,:,:).^2);
     
howmanysamp=numsamp;
for protocolk=1:length(protocolstotest)
    for widthl=1:numel(widths)
        erroratN(protocolk,widthl)=mean(arraytotest(protocolk,widthl,:,howmanysamp));
        errbar(protocolk,widthl)=std(arraytotest(protocolk,widthl,:,howmanysamp));
    end
end
colours={'r','g','b'};

for protocolk=1:length(protocolstotest)
    subplot(1,2,1)
    %CHANGE HERE
%     semilogx(widths,erroratN(protocolk,:))
%     hold on
    shadedErrorBar(widths,erroratN(protocolk,:),errbar(protocolk,:),'lineProps',colours{protocolk});
    h=gca;
    set(h,'xscale','log')
end

%CHANGE HERE
% legend({'Staircase','Bayesian Staircase'})

xlabel('Ratio of step size to \sigma')
ylabel('Error in Mean, %')
ylim([0 inf])


% AGAIN FOR SIGMA
arraytotest=sqrt(errbarsig(:,:,:,:).^2);

for protocolk=1:length(protocolstotest)
    for widthl=1:numel(widths)
        erroratN(protocolk,widthl)=mean(arraytotest(protocolk,widthl,:,howmanysamp));
        errbar(protocolk,widthl)=std(arraytotest(protocolk,widthl,:,howmanysamp));
    end
end
colours={'r','g','b'};

for protocolk=1:length(protocolstotest)
    subplot(1,2,2)
    %CHANGE HERE
%     semilogx(widths,erroratN(protocolk,:))
%     hold on
    shadedErrorBar(widths,erroratN(protocolk,:),errbar(protocolk,:),'lineProps',colours{protocolk});
    h=gca;
    set(h,'xscale','log')
end

%CHANGE HERE
% legend({'Staircase','Bayesian Staircase'})
legend({'Stress Step','Bayesian Stress Step'})

xlabel('Ratio of step size to \sigma')
ylabel('Error in Standard Deviation, %')
ylim([0 inf])

sgtitle(['Convergence after ' num2str(howmanysamp) ' samples'],'FontSize',12)

% % 
% print(['PUBLICATION3 THETASIGConvergence combined error STEP' num2str(nrepeats) ' width ' num2str(distwidth)], '-dpng','-r600')
% savefig(['PUBLICATION3 THETASIGConvergence combined error STEP ' num2str(nrepeats)  ' width  ' num2str(distwidth)])

% print(['PUBLICATION3 THETASIGConvergence combined error STAIR' num2str(nrepeats) ' width ' num2str(distwidth)], '-dpng','-r600')
% savefig(['PUBLICATION3 THETASIGConvergence combined error STAIR ' num2str(nrepeats)  ' width  ' num2str(distwidth)])
