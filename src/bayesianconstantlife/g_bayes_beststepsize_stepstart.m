function [beststress,beststep]=g_bayes_beststepsize_stepstart(theta,sigma,lprior,minstressstep,startingstress)
% Evaluate the next best step size or starting point to test! 
% This is built using a step utility function as described in the paper 
% and in g__STEP_UTILITY
% INPUTS: theta, sigma, log of prior
% OUTPUTS: best stress starting point and step size to test in units of theta

startingstresses=startingstress-minstressstep:2:startingstress;

numstresssteptolookat=1;
if numstresssteptolookat ==1
    stepsizetolookat=minstressstep;
else
    stepsizetolookat=linspace(minstressstep,100,numstresssteptolookat);
end


%% For variable starting point
fun=@(x)g_STEP_UTILITY(minstressstep,x,theta,sigma,lprior);

beststep=minstressstep;

if numel(startingstresses)>40
    %if you have lots of startingstresses:
    % Coarse run through the range of stresses, and see where the utility
    % function actually has values that matter.
    startingstressescoarse=startingstress-minstressstep:10:startingstress;
    z=zeros(1,numel(startingstressescoarse));
    for i=1:numel(startingstressescoarse)
        z(i)=-fun(startingstressescoarse(i));
    end
    
    %In general, it will be a bell shape with 1 or two peaks. 
    culledstresses=startingstressescoarse(z>=(max(z)-0.5*(max(z)-min(z))));
    %DEBUG: look at what's been culled:
    %{
    figure(2)
    plot(startingstressescoarse,z)
    hold on
    plot(culledstresses,z(z>(max(z)-0.5*(max(z)-min(z)))))
    hold off
    %}

    % Now seed based on the culled stresses, 4 points to find the minima. This
    % is fairly efficient, and mostly robust. 
    try
        seedpoints=unique(round(linspace(1,numel(culledstresses),4)));
        opts = optimoptions(@fmincon,'Algorithm','interior-point');
        problem = createOptimProblem('fmincon','objective',...
            fun,'x0',startingstress-0.5*minstressstep,'lb',startingstress-minstressstep,'ub',startingstress,'options',opts);
        ms = MultiStart;
        ms.Display='off';
        tpoints = CustomStartPointSet(culledstresses(seedpoints)');
        [beststress,f] = run(ms,problem,tpoints);
    catch
        z=zeros(1,numel(startingstresses));
        second=zeros(1,numel(startingstresses));
        for i=1:numel(startingstresses)
            [z(i),second(i)]=fun(startingstresses(i));
        end
        z=-1.*z;
        if nnz(isnan(z))>0
            error('Nans in the utility function')
        end
        I=find(z==max(z(:)));
        if nnz(I)>1
            beststress=min(startingstresses(I));
        else
            beststress=startingstresses(I(1));
        end
        f=min(z(I));
    end
else
    %otherwise just do full field since there's a small enough number to do
    %so efficiently
    z=zeros(1,numel(startingstresses));
    second=zeros(1,numel(startingstresses));
    for i=1:numel(startingstresses)
        [z(i),second(i)]=fun(startingstresses(i));
    end
    z=-1.*z;
    if nnz(isnan(z))>0
        error('Nans in the utility function')
    end
    I=find(z==max(z(:)));
    if nnz(I)>1
        beststress=min(startingstresses(I));
    else
        beststress=startingstresses(I(1));
    end

    f=min(z(I));

    %DEBUG
    %{
    figure(2)
    plot(startingstresses,z)
    %}

end

%% For variable step size and starting point 
% At constant life it's always the smallest step size.
%{
fun=@(x)g_STEP_UTILITY(x(1),x(2),theta,sigma,lprior);

NS1=combvec(stepsizetolookat,startingstresses);

z=zeros(1,size(NS1,2));
second=zeros(1,size(NS1,2));
w = waitbar(0);

for i=1:numel(NS1(1,:))
    waitbar(i/numel(NS1(1,:)),w,['Sample number: ' num2str(i),'/',num2str(numel(NS1(1,:)))])
    [z(i),second(i)]=fun([NS1(1,i),NS1(2,i)]);
end
close(w)
z=-1.*z;
if nnz(isnan(z))>0
    error('Nans in the utility function')
end
I=find(z==max(z(:)));
if nnz(I)>1
    [beststep,~]=min(NS1(1,I));
    beststress=min(NS1(2,I));
else
    beststep=NS1(1,I(1));    
    beststress=NS1(2,I(1));
end
%}
%% DEBUG - 
% Plot the utility function
%{
if numstresssteptolookat ==1
    if numel(startingstresses)>40 %if you haven't already done this 
        z=zeros(1,numel(startingstresses));
        second=zeros(1,numel(startingstresses));
        for i=1:numel(startingstresses)
            [z(i),second(i)]=fun(startingstresses(i));
        end

        z=-1.*z;
        if nnz(isnan(z))>0
            error('Nans in the utility function')
        end
    end
    figure(1)
    plot(startingstresses(:),z)
    hold on
    scatter(beststress,max(z));
    ylabel('Utility')
    xlabel('Starting Stress/MPa')
    hold off
    drawnow
else
    gridN=reshape(NS1(1,:),numel(stepsizetolookat), numel(startingstresses));
    gridS=reshape(NS1(2,:),numel(stepsizetolookat), numel(startingstresses));
    gridz=reshape(z,numel(stepsizetolookat), numel(startingstresses));
    h=contourf(gridN,gridS,gridz);
    xlabel('Stress Step /MPa');
    ylabel('Starting Stress /MPa')
    c=colorbar;
    c.Label.String = 'Utility';
    hold on
    I=find(gridz==max(gridz(:)));
    % plot(gridN(I),gridS(I),'ko')
    1;
    hold off
    pause(0.1)
    for i=1:numel(stepsizetolookat)
        figure
        plot(gridz(i,:))
        pause(0.2)
    end
end
%}

% PLOT THE PAPER FIGURE WITH THE VARIATION IN UTILITY FOR SOME NUMBER OF
% SPECIMENS
%{

k=evalin('caller','k');%CMM HACK : copy this line into B_STEPsimulate as well in order for this to work. 
whichsamp=[1,2,3,20,40,70]; %which samples do you want to look at 

if ismember(k,whichsamp)
    z=zeros(1,numel(startingstresses));
    second=zeros(1,numel(startingstresses));
    for i=1:numel(startingstresses)
        [z(i),second(i)]=fun(startingstresses(i));
    end
    z=-1.*z;
    figure(1);
    z=z-min(z);%NORMALISE Z FOR CONVENIENCE
    z=z./max(z);
    z=z-(1+find(k==whichsamp));
    plot(startingstresses,z);
    xlabel('Starting Stress');
    ylabel('Utility Function')
    hold on
    scatter(beststress,max(z),'kx')
    if k==whichsamp(end)
        for i=1:numel(whichsamp)
            legtext{i}=['Sample ',num2str(whichsamp(i))];
        end
        h=gca;
        legend([h.Children(flip(2*(1:numel(whichsamp))))],legtext) %magic....
        set(gca,'Yticklabel', [])
        print(['utilitySTEPfirst',num2str(whichsamp(end))],'-dpng','-r0')
        savefig(['utilitySTEPfirst',num2str(whichsamp(end)),'.fig'])
    end
end

%}
end