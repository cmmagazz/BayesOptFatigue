function totalresults=f_SNresults(totalresults)
% Function to generate SN results
% Current method involves making a stress-step SN curve with full data
% information: the stress and runout of the penultimate test, the stress
% and runout of the last test (failure), and another at that stress with
% the actual cycle to failure. 
% This adds to the totalresults struct with a new field: SN results,
% structure in rows as: penultimate stress, ultimate stress, 
% runout cycle, failure cycle. 
% DEFINITION OF OUTPUT DATA: FAILURESTRESS:
%   1: Failure stress (MPa)
%   2: Number of steps to failure
%   3: Runout Stress (MPa)
%   4: Cycles for Runout Life (log(cycles))
%   5: Cycles to failure (log(cycles))
%   6: Failure Tally (1=fail, 0=runout)
%CMM 2020, RR 2021

%First need to run through totalresults to get the lower bound
samplestested=fieldnames(totalresults);

%If this has been run before, ignore the SN section of totalresults:
idx=find(strcmp('SN',samplestested));
if ~isempty(idx)
    samplestested(idx)=[];
end


%Setup SNresults
SNresults=NaN(numel(samplestested),4);
SNkey=cell(numel(samplestested),2);
for i=1:numel(samplestested)
    %For every sample tested, extract what we care about.
    try
        while ~isfield(totalresults.(samplestested{i}), 'testresults')
            i=i+1;
        end
        testsruni=fieldnames(totalresults.(samplestested{i}).testresults);

        SNresults(i,2)=numel(testsruni);%numbe of steps to failure
        %Pull out the test data
        if numel(testsruni)>1
            %if there are multiple entries
            if strcmp('f',totalresults.(samplestested{i}).testresults.(testsruni{end}).condition)
                %First penultimate test:
                try
                    SNresults(i,3)=totalresults.(samplestested{i}).testresults.(testsruni{end-1}).stress;
                catch
                    SNresults(i,3)=NaN; %if it failed on the first run
                end
                SNresults(i,6)=1;%it failed
                %now the last test
                SNresults(i,1)=totalresults.(samplestested{i}).testresults.(testsruni{end}).stress;
                SNresults(i,4)=log10(totalresults.(samplestested{i}).testresults.(testsruni{end}).cycles);
                SNresults(i,5)=log10(totalresults.(samplestested{i}).testresults.(testsruni{end}).cyclesfail);
            elseif strcmp('r',totalresults.(samplestested{i}).testresults.(testsruni{end}).condition)
                SNresults(i,6)=0;%it ran out
                %If on the other hand there's still only runout data for that:
                SNresults(i,1)=totalresults.(samplestested{i}).testresults.(testsruni{end}).stress;
                %and there's no failure data:
                SNresults(i,4)=log10(totalresults.(samplestested{i}).testresults.(testsruni{end}).cycles);
            end
        else %if there's only 1 test for this anyway:
            %if that was a failure
            if strcmp('f',totalresults.(samplestested{i}).testresults.(testsruni{end}).condition)
                SNresults(i,6)=1;%it failed
                %First penultimate test:
                try
                    SNresults(i,3)=totalresults.(samplestested{i}).testresults.(testsruni{end-1}).stress;
                catch
                    SNresults(i,3)=NaN; %if it failed on the first run
                end
                %now the last test
                SNresults(i,1)=totalresults.(samplestested{i}).testresults.(testsruni{end}).stress;
                SNresults(i,4)=log10(totalresults.(samplestested{i}).testresults.(testsruni{end}).cycles);
                SNresults(i,5)=log10(totalresults.(samplestested{i}).testresults.(testsruni{end}).cyclesfail);
            elseif strcmp('r',totalresults.(samplestested{i}).testresults.(testsruni{end}).condition)
                SNresults(i,6)=0;%it ran out
                %If on the other hand there's still only runout data for that:
                SNresults(i,3)=totalresults.(samplestested{i}).testresults.(testsruni{end}).stress;
                %and there's no failure data:
                SNresults(i,4)=log10(totalresults.(samplestested{i}).testresults.(testsruni{end}).cycles);
            end
        end
        %Pull out a key for reference
        SNkey{i,1}=samplestested{i};
        try
            SNkey{i,2}=totalresults.(samplestested{i}).sampledetails.type;
        catch 
%             warning('No microstructure detected')
        end
    catch
        warning(['No data in sample: ', samplestested{i}])
    end
end
totalresults.SN.results=SNresults;
totalresults.SN.key=SNkey;
end
