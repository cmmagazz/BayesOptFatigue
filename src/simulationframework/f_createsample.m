function meanFS= f_createsample(TestingSet,varargin)
%Monte Carlo method used to create a set of samlpes and 
%"push" this sample away from the mean SN curve based on
%a distribution set: 
%SDistq(1)==1 is gaussian
%SDistq(1)==2 is weibull
%SDistq(1)==3 is log-normal
%For SN models:
%SDistq(2)==1 is Wohler Relationship (default if none given)
%SDistq(2)==2 is Basquin
%SDistq(2)==3 is Bi-linear


% INPUT
%   TestingSet: SDistq, numsamp, basquin, N (life axis), width (MPa)
% OUTPUT
%   TestingSet: meanFS
c=TestingSet.basquin.c;
alpha=TestingSet.basquin.alpha;
beta=TestingSet.basquin.beta;

%% First generate master curves
%Generate the master curve
if numel(TestingSet.SDistq)<2 || TestingSet.SDistq(2)==1
%DEFAULT IS Chris original modified basquin
    meanFS=alpha*(TestingSet.N).^(-beta)+c;
elseif TestingSet.SDistq(2)==2
    %2 = basquin
    meanFS=alpha*(10.^TestingSet.N).^(-beta)+c;

elseif TestingSet.SDistq(2)==3
    %3 = bilinear
    %Let's do bilinear, elbow at 1e6
    [~,Nelbow]=min(abs(TestingSet.N-6));
    %define a line for HCF with C as the limit, and 1e6 at c+100
    NHCF=numel(TestingSet.N)-Nelbow;
    hcfresponse=linspace(c+100,c,NHCF);
    %LCF will be from that point with a gradient of -alpha/5
    lcfresponse=linspace(TestingSet.N(1),6,Nelbow)*(-alpha/5);
    lcfresponse=lcfresponse+(max(hcfresponse)-min(lcfresponse));
    %Generate the master curve
    meanFS=[lcfresponse,hcfresponse];
elseif TestingSet.SDistq(2)==4
    %4 = wohler
    meanFS=-alpha*(TestingSet.N)+c;
else
    error('Model not provided')
end
    
%% Deviate as need be
if TestingSet.SDistq(1)==1
    %Normal random distribution just in stress
    MCF1=normrnd(0,1,[1,TestingSet.numsamp]);
    meanFS=meanFS+MCF1'*TestingSet.width;
elseif TestingSet.SDistq(1)==2
    %weibull random - do elementwise
    MCF1=NaN(TestingSet.numsamp,length(meanFS));
    for i=1:length(meanFS)
        MCF1(:,i)=wblrnd(meanFS(i),TestingSet.width,[1,TestingSet.numsamp]);
    end
    meanFS=MCF1;
elseif TestingSet.SDistq(1)==3
    MCF1=NaN(TestingSet.numsamp,length(meanFS));
    for i=1:length(meanFS)
        MCF1(:,i)=lognrnd(log(meanFS(i)),TestingSet.width,[1,TestingSet.numsamp]);
    end
    meanFS=MCF1;
end

        
%% plot?
if numel(varargin)>0%
    if varargin{1}==1
        semilogx(10.^TestingSet.N,meanFS')
        xlabel('Cycles')
        ylabel('Stress /MPa')
    end
end



end
