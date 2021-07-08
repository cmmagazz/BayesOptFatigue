function samplelevelresults=f_testing(samplelevelresults)
%Function to run fatigue testing for a given sample: easy step process:
%If the sample hasn't been tested before, it takes the runout value and
%asks for a starting stress. If it has been tested before, it applies the
%5% stress step rule (though this is free to change).
%
%Data is structured as follows:
%1) The base details: testdet.runout & *.step are inputted from
%the start: generally these don't change within a sample. *.prevsamp is
%the name of a previous sample that had been tested to that life.
%
%2) The actual test results:To determine the test, you start from the
%front of the alphabet running forwards. This means that the latest
%letter of the alphabet present in testresults is the one to failure /
%latest one. Easy way to do this shown below.

if strcmp(samplelevelresults.testdet.protocol,'stress step')
    %TEMPORARILY PULL THIS OUT, or create it
    try
        testresults=samplelevelresults.testresults;
    catch
        testresults=struct;
    end
    %Setting up the latest test
    s = 'abcdefghijklmnopqrstuvwxyz';
    %First check if it's been tested before:
    i=1;
    while isfield(testresults, s(i))
        i=i+1; %Then want one MORE to make a new test entry
    end
    testresults.(s(i))=struct;
    %The i'th letter is therefore not present and is the one to test. This
    %allows for 26 test.
    %Now calculate the stress that this should be tested at:
    if i>1
        if isfield(testresults.(s(i-1)),'cyclesfail')
            error('Sample already tested to failure')
        end
        %if there was a test WITH THIS SAMPLE, add the step size in MPA
        testresults.(s(i)).stress=testresults.(s(i-1)).stress+samplelevelresults.testdet.step;
        disp(['Next step - increasing stress by ',num2str(samplelevelresults.testdet.step),...
            'MPa']);
        disp(['This means testing to: ',num2str(testresults.(s(i)).stress,4),'MPa']);
    elseif ~isnan(samplelevelresults.testdet.prevsampdet.prevstress) 
        %if there was a previous sample specified, get those results
        % In step testing, you want to start at the same point
        testresults.(s(i)).stress=samplelevelresults.testdet.prevsampdet.prevstress(1);
        disp(['Using first samples first stress: ', num2str(testresults.(s(i)).stress,4),'MPa']);
    else
        %if there's never been a test done:
        if isfield(samplelevelresults.testdet,'startingstress')
            testresults.(s(i)).stress=samplelevelresults.testdet.startingstress;
            disp(['Starting stress set in main script at: ',num2str(testresults.(s(i)).stress,4),'MPa']);
        else
            prompt=('Starting Stress, ideally 70% fatigue limit (MPa): ');
            testresults.(s(i)).stress=input(prompt);
        end
    end
elseif strcmp(samplelevelresults.testdet.protocol,'staircase')
    %TEMPORARILY CHECK THIS OUT, or create it
    if isfield(samplelevelresults,'testresults')
        error('Sample already tested')
    else
        testresults=struct;
    end
    %here we only have one test so 
    s='a';i=1;
    testresults.(s(i))=struct;
    %Now calculate the stress that this should be tested at:
    if ~isnan(samplelevelresults.testdet.prevsampdet.prevstress) 
        %if there was a previous sample specified, get those results
        % In staircase, if it failed go down, if it survived go up
        if samplelevelresults.testdet.prevsampdet.tally==1
            testresults.(s(i)).stress=samplelevelresults.testdet.prevsampdet.prevstress(1)-samplelevelresults.testdet.step;
            disp(['Previous sample failed - decreasing stress by ',num2str(samplelevelresults.testdet.step),...
                'MPa']);
            disp(['This means testing to: ',num2str(testresults.(s(i)).stress,4),'MPa']);
        elseif samplelevelresults.testdet.prevsampdet.tally==0
            testresults.(s(i)).stress=samplelevelresults.testdet.prevsampdet.prevstress(1)+samplelevelresults.testdet.step;
            disp(['Previous sample survived - increase stress by ',num2str(samplelevelresults.testdet.step),...
                'MPa']);
            disp(['This means testing to: ',num2str(testresults.(s(i)).stress,4),'MPa']);
        else
            error('Unable to determine previous sample status')
        end
    else
        %if there's never been a test done:
        if isfield(samplelevelresults.testdet,'startingstress')
            testresults.(s(i)).stress=samplelevelresults.testdet.startingstress;
            disp(['Starting stress set in main script at: ',num2str(testresults.(s(i)).stress,4),'MPa']);
        else
            prompt=('Starting Stress, around your best estimate for mean (MPa): ');
            testresults.(s(i)).stress=input(prompt);
        end
    end
elseif strcmp(samplelevelresults.testdet.protocol,'bayes staircase')
    %TEMPORARILY CHECK THIS OUT, or create it
    if isfield(samplelevelresults,'testresults')
        error('Sample already tested')
    else
        testresults=struct;
    end
    %here we only have one test so 
    s='a';i=1;
    testresults.(s(i))=struct;
    %Now calculate the stress that this should be tested at:
    if ~isnan(samplelevelresults.testdet.prevsampdet.prevstress) 
        %if there was a previous sample specified, get those results
        % In staircase, if it failed go down, if it survived go up
        %create the simulation style failuretally; 

        failurestress=samplelevelresults.testdet.prevsampdet.failurestress;
         % set the upper and lower bounds of possible values of theta here
        if isfield(samplelevelresults.testdet,'theta') 
            theta=samplelevelresults.testdet.theta;
            sigma=samplelevelresults.testdet.sigma;
        else
            disp('Prior space not given, using default values here')
            maxtheta=1000;
            mintheta=100;
            theta=mintheta:1:maxtheta;
            % set upper and lower bound of possible values of sigma here
            maxsigma=400;
            minsigma=1;
            sigma=minsigma:1:maxsigma;
            testresults.(s(i)).stress=B_simulate(failurestress,theta,sigma);
        end
        %calculate the prior if need be
        if isfield(samplelevelresults.testdet,'lprior')
            lprior=samplelevelresults.testdet.lprior;
        else
            lprior=g_calcprior(failurestress,theta,sigma);
        end
        %find the next stress
        testresults.(s(i)).stress=B_simulate(failurestress,theta,sigma,2,lprior);
        disp(['Next sample stress based on Bayesian Staircase at: '...
            ,num2str(testresults.(s(i)).stress,4),'MPa']);
    else
        %if there's never been a test done:
        if isfield(samplelevelresults.testdet,'startingstress')
            testresults.(s(i)).stress=samplelevelresults.testdet.startingstress;
            disp(['Starting stress set in main script at: ',num2str(testresults.(s(i)).stress,4),'MPa']);
        else
            prompt=('Starting Stress, around your best estimate for mean (MPa): ');
            testresults.(s(i)).stress=input(prompt);
        end
    end

elseif strcmp(samplelevelresults.testdet.protocol,'bayes step')
    %TEMPORARILY PULL THIS OUT, or create it
    try
        testresults=samplelevelresults.testresults;
    catch
        testresults=struct;
    end
    s = 'abcdefghijklmnopqrstuvwxyz';
    %First check if it's been tested before:
    i=1;
    while isfield(testresults, s(i))
        i=i+1; %Then want one MORE to make a new test entry
    end
    testresults.(s(i))=struct;
    %The i'th letter is therefore not present and is the one to test. This
    %allows for 26 test.
    %Now calculate the stress that this should be tested at:
    minstressstep=samplelevelresults.testdet.step;

    if i>1
        if isfield(testresults.(s(i-1)),'cyclesfail')
            error('Sample already tested to failure')
        end
        %if there was a test WITH THIS SAMPLE, add the stress step
        testresults.(s(i)).stress=testresults.(s(i-1)).stress+minstressstep;
        disp(['Next step - increasing stress by ',num2str(samplelevelresults.testdet.step),...
            'MPa: ', num2str(testresults.(s(i)).stress,4),'MPa']);
    elseif ~isnan(samplelevelresults.testdet.prevsampdet.prevstress) 
        %if there was a previous sample specified, get those results
        % In staircase, if it failed go down, if it survived go up
        %create the simulation style failuretally; 
        failurestress=samplelevelresults.testdet.prevsampdet.failurestress;
         % set the upper and lower bounds of possible values of theta here
        if isfield(samplelevelresults.testdet,'theta') 
            theta=samplelevelresults.testdet.theta;
            sigma=samplelevelresults.testdet.sigma;
            startingstress=samplelevelresults.testdet.startingstress;
        else
            disp('Prior space not given, using default values here')
            maxtheta=1200;
            mintheta=100;
            theta=mintheta:2:maxtheta;
            % set upper and lower bound of possible values of sigma here
            maxsigma=600;
            minsigma=1;
            sigma=minsigma:1:maxsigma;
            startingstress=samplelevelresults.testdet.startingstress;
        end
        if isfield(samplelevelresults.testdet,'lprior')
            lprior=samplelevelresults.testdet.lprior;
        else
            lprior=g_calcprior(failurestress,theta,sigma);
        end
        testresults.(s(i)).stress=B_STEP_simulate(failurestress,theta,sigma, startingstress,minstressstep,lprior);
        disp(['Next sample stress based on Bayesian Step at: '...
            ,num2str(testresults.(s(i)).stress,4),'MPa']);
    else
        %if there's never been a test done:
        if isfield(samplelevelresults.testdet,'startingstress')
            testresults.(s(i)).stress=samplelevelresults.testdet.startingstress;
            disp(['Starting stress set in main script at: ',num2str(testresults.(s(i)).stress,4),'MPa']);
        else
            prompt=('Starting Stress, ideally 70% fatigue limit (MPa): ');
            testresults.(s(i)).stress=input(prompt);
        end
    end
else
    error('Protocol not supported')
end

%For all methods, save the runout cycles
testresults.(s(i)).cycles=samplelevelresults.testdet.runout;

%% Running the test

prompt=('Runout or fail R/F [R]:');
str=input(prompt,'s');
if ~(strcmpi(str,'R') || strcmpi(str,'F'))  %sanitise input
    str=input(prompt,'s');
end
testresults.(s(i)).condition=str;

%Now deal with the two options:
if strcmpi(testresults.(s(i)).condition,'R')
    %If it was a runout, store the runout cycle condition
    if strcmp(samplelevelresults.testdet.protocol,'step') ...
            || strcmp(samplelevelresults.testdet.protocol,'bayes step')
        disp('RUNOUT - re run the section in order to test again.')
    elseif strcmp(samplelevelresults.testdet.protocol,'staircase') ...
            || strcmp(samplelevelresults.testdet.protocol,'bayes staircase')
        disp('RUNOUT - re run the entire script with a new sample')
    end
elseif strcmpi(testresults.(s(i)).condition,'F')
    %If it failed, store the cycle number and the original runout too
    prompt=('Failure at cycles = ');
    testresults.(s(i)).cyclesfail=input(prompt);
    disp(['FAILURE - at stress ',num2str(testresults.(s(i)).stress),...
        'MPa and runout at ',num2str(testresults.(s(i)).cycles)])
end

%% put results  into the structure
%Finally, order the tests up with the last one at the top:
testresults=orderfields(testresults);
%and put the results back into the struct
samplelevelresults.testresults=testresults;

end