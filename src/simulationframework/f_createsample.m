function meanFS= f_createsample(TestingSet,varargin)
%Monte Carlo method used to create a set of samlpes and 
%"push" this sample away from the mean SN curve based on
%a distribution set: 
%SDistq(1)==1 is gaussian
%SDistq(1)==2 is weibull
%SDistq(1)==3 is log-normal
%For SN models:
%SDistq(2)==1 is Basquin Relationship (default if none given)
%SDistq(2)==2 is Bi-linear

% INPUT
%   TestingSet: SDistq, numsamp, basquin, N (life axis), width (MPa)
% OUTPUT
%   TestingSet: meanFS
c=TestingSet.basquin.c;
alpha=TestingSet.basquin.alpha;
beta=TestingSet.basquin.beta;


if TestingSet.SDistq(1)==1
    MCF1=normrnd(0,1,[1,TestingSet.numsamp]);
    MCF2=normrnd(0,1,[1,TestingSet.numsamp]);
    % have sigma in life be a function of sigma in stress only
    % meanFS=alpha*(TestingSet.N).^(-beta)+c;
    %have sigma in life be linked to sigma in stress by a constant
    %DOUBLE PERTURB PARAMETER: 0 for only perturbing in stress, a value for
    %doing also in life
    doubleperturb=0;
    meanFS=alpha*(TestingSet.N+MCF2'*TestingSet.width.*doubleperturb).^(-beta)+c+MCF1'*TestingSet.width;
elseif TestingSet.SDistq(1)==2
    meanFS=alpha*(TestingSet.N).^(-beta)+c;%generate the master curve
    MCF1=NaN(TestingSet.numsamp,length(meanFS));
    for i=1:length(meanFS)
        MCF1(:,i)=wblrnd(meanFS(i),TestingSet.width,[1,TestingSet.numsamp]);
    end
    meanFS=MCF1;
elseif TestingSet.SDistq(1)==3
    meanFS=alpha*(TestingSet.N).^(-beta)+c;%generate the master curve
    MCF1=NaN(TestingSet.numsamp,length(meanFS));
    for i=1:length(meanFS)
        MCF1(:,i)=lognrnd(log(meanFS(i)),TestingSet.width,[1,TestingSet.numsamp]);
    end
    meanFS=MCF1;
end

%% If you don't want a basquin
if numel(TestingSet.SDistq)>1
    if TestingSet.SDistq(2)==2
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
        if TestingSet.SDistq(1)==1
            MCF1=normrnd(0,1,[1,TestingSet.numsamp]);
            % have sigma in life be a function of sigma in stress only
            %have sigma in life be linked to sigma in stress by a constant
            meanFS=meanFS+MCF1'*TestingSet.width;
        elseif TestingSet.SDistq(1)==2
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
    end
end
        
%% plot?
if numel(varargin)>0%
    if varargin{1}==1
        plot(TestingSet.N,meanFS')
        xlabel('Log(Cycles)')
        ylabel('Stress,MPa')
    end
end



end
