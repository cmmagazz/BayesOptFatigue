function p_HPD(arraytoplot,varargin)
% A function to plot the way that I commonly use for priors. 
% OBLIGATORY INPUTS: 
%       arraytoplot - an array the same size as X and Y to be plotted
%       
% OPTIONAL INPUTS
%       X - an X position array, where X(2,1) is greater than X(1,1)
%           Default is sigma
%       Y - a  Y position array, where Y(1,2) is greater than Y(1,1)
%           Default is theta
%       title - string to give figure a title IN DOUBLE QUOTES
%       limits - array to give limits in colour (e.g. [0 1])
%           Default is mean(arraytoplot) +/- 2 standard deviations
%       units - x/y units, IN DOUBLE QUOTES
%           Default is micron 
%       cunits - C units, placed near colourbar IN DOUBLE QUOTES
%       saveq - 1 if you want to save the figure
%           Default is no
%       resultsdir - default is global resultsdir, but can specify your own
%       IN DOUBLE QUOTES
%       saveasfigq - save as a fig? same as global by default
%
%   Example: XPCcontourf(datastack.H, 'title',"Hardness map", 'cunits',
%   "Hardness /GPa") 
%   This plots the H array in datastack, with the above title and colour
%   units, while assuming the rest as default.
% CMM 2020

try
    theta=evalin('base','theta');
catch
    try
        theta=evalin('base','ResultSet.details.theta');
    catch
        theta=evalin('caller','theta');
    end
end
defaultx=theta{2};
defaulty=theta{1};

defaultCI=0.95;
defaultnewfig=1;

defaultfilename='attime_';
defaultsavefigq=0;

p = inputParser;

validArrayPosNum = @(x) isnumeric(x);
addRequired(p,'arraytoplot',validArrayPosNum);

addOptional(p,'X',defaultx,validArrayPosNum);
addOptional(p,'Y',defaulty,validArrayPosNum);
addOptional(p,'CI',defaultCI,@isscalar);
addOptional(p,'newfig',defaultnewfig,@isnumeric);

addOptional(p,'resultsfilename',defaultfilename,@ischar);
addOptional(p,'savefigq',defaultsavefigq,@isnumeric);

parse(p,arraytoplot,varargin{:});

if p.Results.newfig==1
    figure()
end
maxval=f_HPD(arraytoplot,p.Results.CI);
h=pcolor(p.Results.X,p.Results.Y,exp(arraytoplot));

max2=max(exp(arraytoplot),[],'all');
idmin=exp(arraytoplot)==max2;
[row,col]=find(idmin);

set(h, 'EdgeColor', 'none');
colormap(parula)
hold on

scatter(p.Results.X(col),p.Results.Y(row),'kx')

xlabel('Standard deviation /MPa')
ylabel('Mean Failure Strength /MPa')
caxis([0,maxval])
hold off
c=colorbar;
c.Label.String='Probability';
title('Joint Posterior')

time=clock;
if p.Results.savefigq==1
    print(['Staircaseplot',p.Results.resultsfilename, num2str(time(1:4))],'-dpng','-r0')
    savefig(['Staircaseplot',p.Results.resultsfilename, num2str(time(1:4)),'.fig'])
end
end

