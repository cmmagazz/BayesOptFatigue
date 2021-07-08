function p_SN(failurestress,varargin)
%plot an SN curve with the three relevant points: last runout, failure
%stress, and failure cycle. saves this either in resultsfilename or in
%current folder with a default name
%
%INPUTS: 
%   1) failurestress
%OPTIONAL INPUTS
%   1) resultsfilename, as string DEFAULT: "SNCurve_attimeXXX"
%   2) newfig, 1 for a new figure, 0 for plotting on latest gcf DEFAULT: 1
%   3) savefigq, 1 for save fig, 0 for not DEFAULT: 0
% CMM 2020

defaultfilename='attime_';
defaultnewfig=1;
defaultsavefigq=0;

p = inputParser;

validArrayPosNum = @(x) isnumeric(x);
addRequired(p,'failurestress',validArrayPosNum);

addOptional(p,'resultsfilename',defaultfilename,@ischar);
addOptional(p,'newfig',defaultnewfig,@isnumeric);
addOptional(p,'savefigq',defaultsavefigq,@isnumeric);

parse(p,failurestress,varargin{:});
failurestress=p.Results.failurestress;
%sanitise input
if ~isempty(failurestress)
    idS=isnan(failurestress(:,1)) & isnan(failurestress(:,3));
    failurestress(idS,:)=[];
    %data=f_stripexcessrunouts(data);
end

runout=failurestress(:,4);
failurecycle=failurestress(:,5);

%sanitise the inputs
if nanmean(runout)<100
    runout=10.^runout;
end
if nanmean(failurecycle)<100
    failurecycle=10.^failurecycle;
end


%new fig?
if p.Results.newfig==1
    figure()
end
%plot lower bound (last runout), where it's not nan
p1=semilogx(runout(~isnan(failurestress(:,3))),failurestress(~isnan(failurestress(:,3)),3),'g>');
hold on
%plot upper bound (last test), where it's not nan
p2=semilogx(runout(~isnan(failurestress(:,1))),failurestress(~isnan(failurestress(:,1)),1),'r<');
%plot failurecycle, where not nan
p3=semilogx(failurecycle(~isnan(failurestress(:,1))),failurestress(~isnan(failurestress(:,1)),1),'kx');%cycles to failure
set(0,'DefaultLegendAutoUpdate','off')

if numel(failurecycle(~isnan(failurestress(:,1))))>0 && numel(runout(~isnan(failurestress(:,3))))>0
    %if there are both runouts and failurecycles
    legend([p1,p2,p3], {'Runout', 'Failure', 'Cycle to failure'});
elseif numel(runout(~isnan(failurestress(:,3))))>0
    %if there are no failurecycles, but there are runouts
    legend([p1,p2],{'Runout', 'Failure'});
elseif numel(runout(~isnan(failurestress(:,1))))>0
    %if there are neither runouts nor failurecycles, but there are failures
    %of some sort
    legend([p2],{'Failure'});
end
    
%plot the bars between runouts and failure stresses
if nnz(~isnan(failurestress(:,1)).*~isnan(failurestress(:,3)))>0
    hold on
    for i = 1:size(runout,1)
        plot([runout(i), runout(i)],[failurestress(i,3),failurestress(i,1)],'k')
    end
end
if nnz(~isnan(failurestress(:,1)).*~isnan(failurecycle(:)))>0
    for i=1:numel(failurecycle)
        plot([failurecycle(i), runout(i)],[failurestress(i,1), failurestress(i,1)],'r')
    end
end
xlabel('Cycle Number')
%xlim([0 1e10])
ylabel('Stress /MPa')
title('S/N Diagram')
time=clock;
if p.Results.savefigq==1
    print(['SNCurve',p.Results.resultsfilename, num2str(time(1:4))],'-dpng','-r0')
end
end