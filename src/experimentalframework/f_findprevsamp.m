function prevsampdet=f_findprevsamp(totalresults,testdet,varargin)
%Simple function to just find the previous sample at that life
runoutvalue=testdet.runout; %runout value for this sample
stressesdone=[];
%find all samples tested
samplestested=fieldnames(totalresults);

%If SN results has been run before, ignore the SN section of totalresults:
idx=find(strcmp('SN',samplestested));
if ~isempty(idx)
    samplestested(idx)=[];
end
%Now remove all samples that had been done in a different protocol
idx=NaN(numel(samplestested),1);
for i=1:numel(samplestested)
    try
        if ~strcmp(testdet.protocol,totalresults.(samplestested{i}).testdet.protocol)
            idx(i)=i;
        end
    catch
        idx(i)=i;%if it doesn't even have the protocol saved
    end
end
idx(isnan(idx))=[];
if ~isempty(idx)
    samplestested(idx)=[];
end

if numel(samplestested)==0
    warning('No samples detected OR not valid protocol selected')
end

%Go through each previously done sample and see they are the same protocol.
%If they are, use that information to get the previous samples history, so
%that you can make an informed decision on the next sample.
%loop through in order to get the stress at which it was tested IF the
%runout values were the same
if strcmp(testdet.protocol,'staircase')
    if numel(samplestested)> 1
        %need to find the last sample:
        lastsample=samplestested{end-1};
        testsruni=fieldnames(totalresults.(lastsample).testresults);
        if totalresults.(lastsample).testdet.runout==runoutvalue
            stressesdone=totalresults.(lastsample).testresults.(testsruni{end}).stress;%get the last stress
            if strcmp(totalresults.(lastsample).testresults.(testsruni{end}).condition,'f')
                tallyi=1; %if it failed, set to 1
            else
                tallyi=0;
            end
        end
    end
    %if there were any, find the mean and stdev.
    if numel(samplestested)>1
        prevsampdet.prevstress=stressesdone;
        prevsampdet.tally=tallyi;
    else
        prevsampdet.prevstress=NaN;
        prevsampdet.tally=NaN;
    end
elseif strcmp(testdet.protocol,'bayes staircase')
    %need to find the last sample:
    stressesdone=[];
    failuretally=[];
    
    if numel(samplestested)> 1
        failurestress=NaN(numel(samplestested)-1,6);
        for i=1:numel(samplestested)-1     %pick all but the current sample
            testsruni=fieldnames(totalresults.(samplestested{i}).testresults);% take the last test
            if totalresults.(samplestested{i}).testdet.runout==runoutvalue %if it was at that runout value
                stressesdone=[stressesdone, totalresults.(samplestested{i}).testresults.(testsruni{end}).stress];%get the last stress
                if strcmp(totalresults.(samplestested{i}).testresults.(testsruni{end}).condition,'f')
                    tallyi=1; %if it failed, set to 1
                    failurestress(i,:)=[stressesdone(end),NaN,NaN,totalresults.(samplestested{i}).testdet.runout,totalresults.(samplestested{i}).testresults.(testsruni{end}).cyclesfail,1]; %if it failed, set to 1
                else
                    tallyi=0;
                    failurestress(i,:)=[NaN,NaN,stressesdone(end),totalresults.(samplestested{i}).testdet.runout,NaN,0]; %if it failed, set to 1
                end
                failuretally=[failuretally, tallyi];
            end
        end
    end
    if numel(samplestested)>1
        prevsampdet.prevstress=stressesdone;
        prevsampdet.tally=failuretally;
        prevsampdet.failurestress=failurestress;
    else
        prevsampdet.prevstress=NaN;
        prevsampdet.tally=NaN;
        prevsampdet.failurestress=NaN;
    end
elseif strcmp(testdet.protocol,'stress step')
    if numel(samplestested)> 1
        %need to find the last sample:
        lastsample=samplestested{end-1};
        testsruni=fieldnames(totalresults.(lastsample).testresults);
        if totalresults.(lastsample).testdet.runout==runoutvalue
            stressesdone=totalresults.(lastsample).testresults.(testsruni{1}).stress;%get the last stress
            if strcmp(totalresults.(lastsample).testresults.(testsruni{1}).condition,'f')
                tallyi=1; %if it failed, set to 1
            else
                tallyi=0;
            end
        end
    end
    %if there were any, find the mean and stdev.
    if numel(stressesdone)>0
        prevsampdet.prevstress=stressesdone;
        prevsampdet.tally=tallyi;
    else
        prevsampdet.prevstress=NaN;
        prevsampdet.tally=NaN;
    end

elseif strcmp(testdet.protocol,'bayes step')
    %need to find the last sample details:
    stressesdone=[];
    penultimatestress=[];
    failuretally=[];
    if numel(samplestested)> 1
        failurestress=NaN(numel(samplestested)-1,6);
        for i=1:numel(samplestested)-1     %pick all but the current sample
            testsruni=fieldnames(totalresults.(samplestested{i}).testresults);% take the last test
            if totalresults.(samplestested{i}).testdet.runout==runoutvalue %if it was at that runout value
                stressesdone=[stressesdone, totalresults.(samplestested{i}).testresults.(testsruni{end}).stress];%get the last stress
                try
                    penultimatestress=[penultimatestress,totalresults.(samplestested{i}).testresults.(testsruni{end-1}).stress];
                catch
                    penultimatestress=[penultimatestress,NaN];
                end
                if strcmp(totalresults.(samplestested{i}).testresults.(testsruni{end}).condition,'f')
                    failurestress(i,:)=[stressesdone(end),numel(testsruni),penultimatestress(end),totalresults.(samplestested{i}).testdet.runout,totalresults.(samplestested{i}).testresults.(testsruni{end}).cyclesfail,1]; %if it failed, set to 1
                    failuretally=[stressesdone, 1];
                else
                    failurestress(i,:)=[NaN,numel(testsruni),stressesdone(end),totalresults.(samplestested{i}).testdet.runout,NaN,0]; %if it failed, set to 1
                    failuretally=[stressesdone, 0];
                end
            end
        end
    end
    if numel(samplestested)>1
        prevsampdet.prevstress=stressesdone;
        prevsampdet.failurestress=failurestress;
        prevsampdet.tally=failuretally;
    else
        prevsampdet.prevstress=NaN;
        prevsampdet.failurestress=NaN;
        prevsampdet.tally=NaN;
    end

else
    error('Protocol not supported')
end
end