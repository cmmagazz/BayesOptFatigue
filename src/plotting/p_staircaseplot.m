function p_staircaseplot(failurestress,varargin)
%plot an staircase plot with a failure tally input. 
% saves this either in resultsfilename or in
% current folder with a default name
%
%INPUTS: 
%   EITHER) FAILURESTRESS as the 6 column definition
%   OR) FAILURETALLY as the 2 column definition
%Optional:
%   'resultsfilename': String for filename
%   'savefigq': 1 save, 0 dont save

defaultfilename='attime_';
defaultsavefigq=0;

p = inputParser;

validArrayPosNum = @(x) isnumeric(x);
addRequired(p,'failurestress',validArrayPosNum);
addOptional(p,'resultsfilename',defaultfilename,@ischar);
addOptional(p,'savefigq',defaultsavefigq,@isnumeric);

parse(p,failurestress,varargin{:});

%Figure out the top point for all of failurestress regardless of step or
%staircase
topdata=NaN(size(failurestress,1),2);

if size(failurestress,2)==2
    topdata=failurestress;
elseif size(failurestress,2)==6
    topdata=[max([failurestress(:,1),failurestress(:,3)],[],2),failurestress(:,6)];
end
 
figure();
sampleN=1:size(topdata,1);
idrunout=topdata(:,2)==0;
idfail=topdata(:,2)==1;
failuretallyF=topdata(:,1);
failuretallyS=failurestress(idrunout,3);

failuretallyF=failuretallyF(idfail);
%Draw the top points
scatter(sampleN(idfail),failuretallyF,'x','MarkerEdgeColor',[0 0 0])
hold on

%Now draw the runouts: first consider the runouts that were there
%registered directly in failurestress, linked to a sample number
if size(failurestress,2)==2
    failuretallyS=[sampleN(idrunout),failurestress(idrunout,1)];
elseif size(failurestress,2)==6
    failuretallyS=[sampleN(idrunout)',failurestress(idrunout,3)];
    %Now add in the steps from any possible staircase test
    for i=1:numel(sampleN)
        numsteps=failurestress(i,2);
        if numsteps>1
            stepsize=failurestress(i,1)-failurestress(i,3);
            runoutstressesforsampleI=failurestress(i,1)-([2:numsteps]-1)*stepsize;
            failuretallyS=[failuretallyS;repmat(sampleN(i),(numsteps-1),1),runoutstressesforsampleI'];
        end
    end
end
%now add these in
scatter(failuretallyS(:,1),failuretallyS(:,2),'o','MarkerEdgeColor',[0 0 0])
legend({'Failure','Runout'},'AutoUpdate','off')

plot(sampleN,topdata(:,1),'k')
if ~strcmp(p.Results.resultsfilename,'attime_')
    title(p.Results.resultsfilename)
else
    title('Staircase Diagram')
end

xlabel('Sample number')
ylabel('Test Stress /MPa')
if max(sampleN)<10
    xticks(1:max(sampleN))
end
hold off
time=clock;

if p.Results.savefigq==1
    print(['Staircaseplot',p.Results.resultsfilename, num2str(time(1:4))],'-dpng','-r0')
    savefig(['Staircaseplot',p.Results.resultsfilename, num2str(time(1:4)),'.fig'])
end
end