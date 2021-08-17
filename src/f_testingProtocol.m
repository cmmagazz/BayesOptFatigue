function [ResultSetraw]=f_testingProtocol(TestingSet,ResultSet)
%Take the testing set data and perform the desired protocol to it and
%save the results back into result set. Each protocol is a different
%way of testing the previousley generated sample SN curves and following the stress
%and evaluating if it fails
%
%PROTOCOLS AVAILABLE:
%       1: Standard step testing
%       2: Staircase
%       3: Probit
%       4: Life testing
%       5: Bayesian Staircase
%       6: Bayesian Step
%INPUTS:
%   TestingSet: in one go
%   ResultSet: in one go
%
%OUTPUTS:
%   ResultSetraw: ONLY THE RAW branch of the struct - but in reality the only things changed are:
%       failurestress is array of fatigue results that is relevant.
%       shannon if necessary
%       the prior if necessary

% DEFINITION OF OUTPUT DATA: FAILURESTRESS AS:
%   1: Failure stress (MPa)
%   2: Number of steps to failure
%   3: Runout Stress (MPa)
%   4: Cycles for Runout Life (log(cycles))
%   5: Cycles to failure (log(cycles))
%   6: Failure Tally (1=fail, 0=runout)
numsamp = TestingSet.details.numsamp;
try
    startingstress=ResultSet.details.startingstress;
catch
    ResultSet.details=ResultSet;
    fn=fieldnames(ResultSet);
    ResultSet = rmfield(ResultSet,fn(~strcmp(fn,'details')));
    startingstress=ResultSet.details.startingstress;
end

stepsize=ResultSet.details.step.stepsize;


failurestress=NaN(numsamp,6);

if strcmp(ResultSet.details.protocol,'stress step')
    for k=1:numsamp
        %the classic stress step: from a given starting point, increase
        %until it breaks
        [failurestress(k,1),failurestress(k,2),failurestress(k,3),failurestress(k,5)]=f_testingDAM(TestingSet.details,TestingSet.meanFS(k,:),ResultSet.details);
        failurestress(k,6)=1;
        failurestress(k,4)=ResultSet.details.runout;
        
        if isfield(ResultSet, 'plotq')
            if ResultSet.plotq==1
                figure(2)
                p_SN(failurestress,'newfig',0)
                drawnow
            end
        end
    end
elseif strcmp(ResultSet.details.protocol,'staircase')
    %staircase method described in early literature - if the sample failes
    %go down a stress level for the next sample and if it survives go up a
    %stress level for the next sample.
    teststresses=startingstress;
    stepsizeabs=ResultSet.details.step.stepsize;
    for k=1:numsamp
        FSsample=TestingSet.meanFS(k,:); %pick out the sample
        %element where life=runout number
        [~,idN]=min(abs(TestingSet.details.N(:)-ResultSet.details.runout));
        condition=floor(teststresses/FSsample(idN));
        condition(condition>1)=1;
        condition(condition<0)=0;
        failurestress(k,6)=condition;
        
        if condition==1 %if it failed
            failurestress(k,1)=teststresses;
            %additionally, find the cycle at which it failed
            [~,idS]=min(abs(FSsample(:)-teststresses));
            failurestress(k,5)=TestingSet.details.N(idS);
            teststresses=teststresses-stepsizeabs;
        elseif condition==0
            %if it didn't fail
            failurestress(k,3)=teststresses;
            teststresses=teststresses+stepsizeabs;
        end
        failurestress(k,4)=ResultSet.details.runout;
        
        if isfield(ResultSet, 'plotq')
            if ResultSet.plotq==1
                figure(2)
                p_SN(failurestress,'newfig',0)
                drawnow
            end
        end
    end
    
elseif strcmp(ResultSet.details.protocol,'probit')
    %Probit method described in early literature - a range of stresses
    %is defined and the samples are allocated to a bin corresponding
    %to one of the stresses.
    %disp(['For Probit method, using number of step sizes as number of stress steps, i.e. ', num2str(ResultSet.details.step.numsizes)])

    %using the ASTM definition of around 5 stresses:
    numstresses=5;
    stresses=startingstress:stepsize:startingstress+stepsize*numstresses;
    [~,idN]=min(abs(TestingSet.details.N(:)-ResultSet.details.runout));
    %divide the samples into buckets for each stress
    bucketsize=floor(TestingSet.details.numsamp/size(stresses,2));
    if bucketsize>1 && max(stresses)>(mean(TestingSet.meanFS(:,idN))-std(TestingSet.meanFS(:,idN)))
        failuretallytemp=zeros(size(stresses,2),bucketsize);
        for j=1:size(stresses,2)
            samplestotest=TestingSet.meanFS((j-1)*bucketsize+1:j*bucketsize,:);
            for k=1:size(samplestotest,1)%built like the others, but simpler
                FSsample=samplestotest(k,:); %pick out the sample
                
                %element where life=runout number
                [~,idN]=min(abs(TestingSet.details.N(:)-ResultSet.details.runout));
                
                condition=floor(stresses(j)/FSsample(idN));
                condition(condition>1)=1;
                failuretallytemp(j,k)= condition;
                if condition==1    %additionally, find the cycle at which they all failed
                    [~,idS]=min(abs(FSsample(:)-stresses(j)));
                    failurestress(:,k+(j-1)*size(samplestotest,1))=TestingSet.details.N(idS);
                end
            end
        end
        failuretally=[repelem(stresses',bucketsize),reshape(failuretallytemp',[],1 )]; %arrange in the same manner as others
        failurestress(failuretally(:,2)==1,1)=failuretally(failuretally(:,2)==1,1);
        failurestress(failuretally(:,2)==0,3)=failuretally(failuretally(:,2)==0,1);
        failurestress(1:numel(stresses)*bucketsize,6)=failuretally(:,2);
        failurestress(:,4)=ResultSet.details.runout*ones(numsamp,1);
    else
        warning('Probit set up to get insufficient data.')
    end
    
elseif strcmp(ResultSet.details.protocol,'life')
    %The standard method normally used - life testing.
    %Each sample is tested at constant stress for a range of stresses until
    %failure
    for k=1:size(TestingSet.meanFS,1)%built like the others, but simpler
        FSsample=TestingSet.meanFS(k,:); %now want to find the life at that stress
        [~,idN]=min(abs(FSsample(:)-startingstress));
        failurestress(k,5)=TestingSet.details.N(idN);
        if failurestress(k,5)<ResultSet.details.runout %if this is less than runout, then it failed
            failurestress(k,6)=1;
            failurestress(k,1)=startingstress;
            [~,idS]=min(abs(FSsample(:)-startingstress));
            failurestress(k,5)=TestingSet.details.N(idS);
        else
            failurestress(k,6)=0;
            failurestress(k,3)=startingstress;
        end
        failurestress(k,4)=ResultSet.details.runout;
        
        if isfield(ResultSet, 'plotq')
            if ResultSet.plotq==1
                figure(2)
                p_SN(failurestress,'newfig',0)
                drawnow
            end
        end
    end

elseif strcmp(ResultSet.details.protocol,'bayes staircase')||strcmp(ResultSet.details.protocol,'bayes staircase flat prior')...
        ||strcmp(ResultSet.details.protocol,'bayes staircase 1/x prior')
   
   theta=ResultSet.details.theta;

    if strcmp(ResultSet.details.protocol,'bayes staircase 1/x prior')
        lpriorq=1;%do you want to initialize the prior (1 for 1/x, 2 for flat)
    elseif strcmp(ResultSet.details.protocol,'bayes staircase flat prior')||strcmp(ResultSet.details.protocol,'bayes staircase')
        lpriorq=2;
    end
    lprior=g_calcprior([],theta,lpriorq); %initialise the prior
    z = waitbar(0);
    shannon=NaN(numsamp,1);
    ResultSet.raw.lprior=lprior;
    for k=1:numsamp
        waitbar(k/numsamp,z,['Sample number: ' num2str(k),'/',num2str(numsamp)])
        teststresses=g_bayesbeststress(ResultSet);
        
        FSsample=TestingSet.meanFS(k,:); %pick out the sample
        
        %element where life=runout number
        [~,idN]=min(abs(TestingSet.details.N(:)-ResultSet.details.runout));
        condition=floor(teststresses/FSsample(idN));
        if condition>=1 %if it failed
            failurestress(k,6)=1;
            failurestress(k,1)=teststresses;
            %additionally, find the cycle at which it failed
            [~,idS]=min(abs(FSsample(:)-teststresses));
            failurestress(k,5)=TestingSet.details.N(idS);
        elseif condition==0
            failurestress(k,6)=0;
            failurestress(k,3)=teststresses;
        end
        failurestress(k,4)=ResultSet.details.runout;
        
        ResultSet.raw.failurestress=failurestress;
        [lprior,~,shannon(k)]=g_calcprior(ResultSet,[],lprior,lpriorq);
        ResultSet.raw.lprior=lprior;
        if isfield(ResultSet, 'plotq')
            if ResultSet.plotq==1
                figure(2)
                subplot(1,2,1)
                p_SN(failurestress,'newfig',0)
                subplot(1,2,2)
                p_HPD(lprior,'newfig',0)
                drawnow
            end
        end
    end

    close(z)
    [lprior,~,shannon(k+1)]=g_calcprior(failurestress,theta,lpriorq,lprior);
    ResultSetraw.shannon=shannon;
    ResultSetraw.lprior=lprior;
elseif strcmp(ResultSet.details.protocol,'bayes step')
    theta=ResultSet.details.theta;

    minstressstep=ResultSet.details.step.stepsize;
    origstartstress=ResultSet.details.startingstress;

    [newstartstress,beststep,lprior,shannon]=B_STEP_simulate([],theta,origstartstress,minstressstep);
    z = waitbar(0);

    for k=1:numsamp
        waitbar(k/numsamp,z,['Sample number: ' num2str(k),'/',num2str(numsamp)])
        ResultSet.details.startingstress=newstartstress;
        ResultSet.details.step.stepsize=beststep;
        ResultSetraw.shannon(k)=shannon;  
        [failurestress(k,1),failurestress(k,2),failurestress(k,3),failurestress(k,5)]=f_testingDAM(TestingSet.details,TestingSet.meanFS(k,:),ResultSet.details);
        failurestress(k,6)=1;
        [newstartstress,beststep,lprior,shannon]=B_STEP_simulate(failurestress,theta,origstartstress,minstressstep,lprior);
        failurestress(k,4)=ResultSet.details.runout;
        if isfield(ResultSet, 'plotq')
            if ResultSet.plotq==1
                figure(2)
                p_SN(failurestress,'newfig',0)
                drawnow
            end
        end
    end
    [lprior,~,shannon,~]=g_calcprior(failurestress(end,:),theta,lprior);
    close(z)

    ResultSetraw.shannon(k+1)=shannon;
    ResultSetraw.lprior=lprior;
end
% Save to results
ResultSetraw.failurestress=failurestress;
end
