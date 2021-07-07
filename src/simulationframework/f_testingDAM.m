function [failurestress,failurestep,penultimatestress,failurecycle]=f_testingDAM(TestingSetdetails,SNsample, ResultSetdetails)
%Function to induce, if needed, damage into the sample, should the same
%sample be tested at multiple stresses.
%
% MODELS AVAILABLE:
%       1: No Damage
%       2: Miner's Rule
%
% INPUTS:
%   SNsample: sn curve of 1 sample (the one being tested)
%   TestingSet.details: in one go
%   ResultSet.details: in one go
%
% OUTPUTS:
% failurestress: the stress that it failed at
% failurestep: the step that it failed at


%deal with inputs:
DAMq=TestingSetdetails.DAMq;
N=TestingSetdetails.N;
startingstress=ResultSetdetails.startingstress;
runout=ResultSetdetails.runout;
stepsize=ResultSetdetails.step.stepsize;
maxSN=max(SNsample(:));

%run the script
numsteps=10.*ceil(((maxSN+100)-startingstress)/stepsize);%maximum number of steps needed to reach top of sn curve.
try
    teststresses = NaN(numsteps,1);
catch %if there are no steps at all, then just do 5 
    teststresses = NaN(100,1);
end
teststresses(1)=startingstress;
for i=2:numsteps
    if teststresses(i-1)+stepsize>maxSN+100
        numsteps=i-1; %if you've reached the top of the sn curve, stop.
        break
    end
    teststresses(i)=teststresses(i-1)+stepsize;
end

idNaN=isnan(teststresses);
teststresses(idNaN)=[];
%element where life=runout number
[~,idN]=min(abs(N(:)-runout));

if DAMq==1
    failurestress=interp1(teststresses,teststresses,SNsample(idN),'next');
    if isnan(failurestress)
        failurestress=teststresses(1);
    end
    failurestep=find(teststresses==failurestress);
    [~,idS]=min(abs(SNsample(:)-failurestress));
    failurecycle=TestingSetdetails.N(idS);

elseif DAMq==2
    %case 2: miner's rule. here, we add a damage factor based on how much of a
    %proportion of its life it survived at the previous tests.
    %need to find number of cycles to failure at a given stress:
    condition=0;
    lifeatthatstress=NaN(numsteps,1);
    for i=1:numsteps
        [~,idS]=min(abs(SNsample(:)-teststresses(i)));
        lifeatthatstress(i)=N(idS);
        if lifeatthatstress(i)==max(N) %if at the end of the SN curve, assume fatigue limit
            lifeatthatstress(i)=inf;
        end
        condition=(10^runout)/(10^lifeatthatstress(i))+condition;
        
        if condition>1
            failurestep=i; %find the step that this occured in
            failurestress=teststresses(failurestep); %find the stress of that step
            
            %additionally, find the cycle at which it failed
            %use miner's to find how much of a life debit it would have
            if i==1
                failurecycle=lifeatthatstress(i);
            elseif i>1
                failurecycle=runout*(1-runout/nansum(lifeatthatstress(lifeatthatstress<inf)));
            end
            break
        end
    end
end

if failurestep>1
    penultimatestress=teststresses(failurestep-1);
else
    penultimatestress=NaN;
end

end