function [c,h]=p_contourHPD(arraytoplot,varargin)
% A function to plot the contour and create line handles. 
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

try
    defaultx=evalin('base','sigma');
    defaulty=evalin('base','theta');
catch
    try
        defaultx=evalin('base','ResultSet.details.sigma');
        defaulty=evalin('base','ResultSet.details.theta');
    catch
        defaultx=evalin('caller','sigma');
        defaulty=evalin('caller','theta');
    end
end
defaultCI=[0.5,0.9];
defaultnewfig=1;

p = inputParser;

validArrayPosNum = @(x) isnumeric(x);
addRequired(p,'arraytoplot',validArrayPosNum);

addOptional(p,'X',defaultx,validArrayPosNum);
addOptional(p,'Y',defaulty,validArrayPosNum);
addOptional(p,'CI',defaultCI,@isscalar);
addOptional(p,'newfig',defaultnewfig,@isnumeric);


parse(p,arraytoplot,varargin{:});

if p.Results.newfig==1
    figure()
end
maxval=NaN(length(p.Results.CI),1);
for i=1:length(p.Results.CI)
    maxval(i)=f_HPD(arraytoplot,p.Results.CI(i));
end
[c,h]=contour(p.Results.X,p.Results.Y,exp(arraytoplot),flipud(maxval));

max2=max(exp(arraytoplot),[],'all');
idmin=exp(arraytoplot)==max2;
[row,col]=find(idmin);
hold on
scatter(p.Results.X(col),p.Results.Y(row),'kx')

xlabel('\sigma')
ylabel('\theta')
hold off
end

