%% Divide totalresults into 4 sets based on protocol
staircaseresults=totalresults;
samplestested=fieldnames(staircaseresults);

%If SN results has been run before, ignore the SN section of totalresults:
idx=find(strcmp('SN',samplestested));
if ~isempty(idx)
    samplestested(idx)=[];
end
%Now remove all samples that had been done in a different protocol
idx=NaN(numel(samplestested),1);
for i=1:numel(samplestested)
    if ~strcmp('staircase',staircaseresults.(samplestested{i}).testdet.protocol)
        staircaseresults = rmfield(staircaseresults,samplestested{i});
    end
end

bayesresults=totalresults;
samplestested=fieldnames(bayesresults);
%If SN results has been run before, ignore the SN section of totalresults:
idx=find(strcmp('SN',samplestested));
if ~isempty(idx)
    samplestested(idx)=[];
end
%Now remove all samples that had been done in a different protocol
idx=NaN(numel(samplestested),1);
for i=1:numel(samplestested)
    if ~strcmp('bayes staircase',bayesresults.(samplestested{i}).testdet.protocol)
        bayesresults = rmfield(bayesresults,samplestested{i});
    end
end

stepresults=totalresults;
samplestested=fieldnames(stepresults);
%If SN results has been run before, ignore the SN section of totalresults:
idx=find(strcmp('SN',samplestested));
if ~isempty(idx)
    samplestested(idx)=[];
end
%Now remove all samples that had been done in a different protocol
idx=NaN(numel(samplestested),1);
for i=1:numel(samplestested)
    if ~strcmp('step',stepresults.(samplestested{i}).testdet.protocol)
        stepresults = rmfield(stepresults,samplestested{i});
    end
end


bayesstepresults=totalresults;
samplestested=fieldnames(bayesstepresults);
%If SN results has been run before, ignore the SN section of totalresults:
idx=find(strcmp('SN',samplestested));
if ~isempty(idx)
    samplestested(idx)=[];
end
%Now remove all samples that had been done in a different protocol
idx=NaN(numel(samplestested),1);
for i=1:numel(samplestested)
    if ~strcmp('bayes step',bayesstepresults.(samplestested{i}).testdet.protocol)
        bayesstepresults = rmfield(bayesstepresults,samplestested{i});
    end
end

%generate the sn data
staircaseresults=f_SNresults(staircaseresults);

stepresults=f_SNresults(stepresults);
bayesresults=f_SNresults(bayesresults);
bayesstepresults=f_SNresults(bayesstepresults);

%% convergence plot
%construct failuretallies

%calculate priors after every number of samples
maxtheta=900;
mintheta=10;
theta=mintheta:1:maxtheta;
% set upper and lower bound of possible values of sigma here
maxsigma=300;
minsigma=1;
sigma=minsigma:1:maxsigma;
prior=cell(1,4);
progthet=zeros(30,12);
progsig=zeros(30,12);

[prior{1},~,~,~,progthet(:,1:3),progsig(:,1:3)]=g_calcprior(stepresults.SN.results,{theta,sigma});
[prior{2},~,~,~,progthet(:,4:6),progsig(:,4:6)]=g_calcprior(bayesresults.SN.results,{theta,sigma});
[prior{3},~,~,~,progthet(:,7:9),progsig(:,7:9)]=g_calcprior(staircaseresults.SN.results,{theta,sigma});
[prior{4},~,~,~,progthet(:,10:12),progsig(:,10:12)]=g_calcprior(bayesstepresults.SN.results,{theta,sigma});

%%
%Plot each of the lines 
subplot(1,2,1)
shadedErrorBar(1:30,progthet(:,7),[abs(progthet(:,8)-progthet(:,7)),abs(progthet(:,9)-progthet(:,7))],'lineProps','b')
hold on
shadedErrorBar(1:30,progthet(:,4),[abs(progthet(:,5)-progthet(:,4)),abs(progthet(:,6)-progthet(:,4))],'lineProps','g')
shadedErrorBar(1:30,progthet(:,1),[abs(progthet(:,2)-progthet(:,1)),abs(progthet(:,3)-progthet(:,1))],'lineProps','r')
shadedErrorBar(1:30,progthet(:,10),[abs(progthet(:,11)-progthet(:,10)),abs(progthet(:,12)-progthet(:,10))],'lineProps','m')
legend({'Staircase','Bayesian Staircase','Stress Step','Bayes Step'})
xlabel('Number of Samples')
ylabel('Estimate of Mean')

subplot(1,2,2)
shadedErrorBar(1:30,progsig(:,7),[abs(progsig(:,8)-progsig(:,7)),abs(progsig(:,9)-progsig(:,7))],'lineProps','b')
hold on
shadedErrorBar(1:30,progsig(:,4),[abs(progsig(:,5)-progsig(:,4)),abs(progsig(:,6)-progsig(:,4))],'lineProps','g')
shadedErrorBar(1:30,progsig(:,1),[abs(progsig(:,2)-progsig(:,1)),abs(progsig(:,3)-progsig(:,1))],'lineProps','r')
shadedErrorBar(1:30,progsig(:,10),[abs(progsig(:,11)-progsig(:,10)),abs(progsig(:,12)-progsig(:,10))],'lineProps','m')
legend({'Staircase','Bayesian Staircase','Stress Step','Bayes Step'})
xlabel('Number of Samples')
ylabel('Estimate of Standard Deviation')


print(['Experimental_parameter_convergence_v2'], '-dpng','-r1500')