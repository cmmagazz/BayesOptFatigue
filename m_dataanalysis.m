%% settings
filename='templatefatiguedata1.xlsx';
IDcol='ID';
stresscol='σa / Mpa';
Nfcol='Nf';
Numcycles=5e4;

%% read the data and put into correct formatt
totalresults=struct();
data=readtable(filename,'VariableNamingRule','preserve');
failrows=data.(Nfcol)<Numcycles;
runoutrows=data.(Nfcol)>Numcycles;
ResultSet=struct();
ResultSet.details.theta{1}=0:0.1:150;
ResultSet.details.theta{2}=1:1:1000;
ResultSet.details.dist='2pwbl';
%put the table into the correct format to be used by the other functions
ResultSet.raw.failurestress=NaN([length(failrows),6]);
ResultSet.raw.failurestress(failrows,1)=data.(stresscol)(failrows);
ResultSet.raw.failurestress(runoutrows,3)=data.(stresscol)(runoutrows);
ResultSet.raw.failurestress(:,4)=ones([length(failrows),1]).*log10(Numcycles);
ResultSet.raw.failurestress(failrows,5)=data.(Nfcol)(failrows);
ResultSet.raw.failurestress(failrows,6)=1;
ResultSet.raw.failurestress(runoutrows,6)=0;
%% Prior plots
figure()
lprior=g_calcprior(ResultSet,2);
% set upper and lower bound of possible values of sigma here
p_HPD(lprior)

%% Joint probability plot
figure()
lpriorstack{1}=lprior;%cell array of the joint distributions we are plotting
leg='results';%legend label for the series
xl=[0,200];%xlimits for the plot
yl=[100,150];%ylimits for the plot
p_priorcomparison(lpriorstack,ResultSet,leg,xl,yl)

%% marginal probability CDF plot
figure()
P=NaN(1,length(ResultSet.details.theta{1}));
for i=1:length(ResultSet.details.theta{1})
    temp=g_calcprobCDF(ResultSet.details.theta{1}(i),ResultSet.details.theta,ResultSet.details.dist);
    P(i)=sum(temp.*exp(lprior),'all','omitnan');
end
pPDF=diff(P)./diff(ResultSet.details.theta{1});
pPDF=[pPDF,0];
plot(ResultSet.details.theta{1},pPDF)
xlim([100,150])