function beststress=g_bayesbeststress(ResultSet)
% Evaluate the next best stress to test! This is built using a utility
% function as described in the paper and in g_UTILITY
% INPUTS: theta, sigma, log of prior
% OUTPUTS: best stress to test in units of theta
theta=ResultSet.details.theta;
lprior=ResultSet.raw.lprior;

fun=@(stress)g_UTILITY(stress,ResultSet);
minmu=min(theta{1});
maxmu=max(theta{1});
rangmu=maxmu-minmu;

%% FIRST PASS
% Coarse run through the range of stresses, and see where the utility
% function actually has values that matter. 
stresstolookat=minmu:30:maxmu; %COARSE
z=stresstolookat;
for i=1:numel(stresstolookat)
    z(i)=-fun(stresstolookat(i));
end
%In general, it will be a bell shape with 1 or two peaks. 
%Now we can cull all the stresses that lead to a utility that's less than
%half the range, which will include all the flat bits. 
culledstresses=stresstolookat(z>(max(z)-0.5*(max(z)-min(z))));
%DEBUG: look at what's been culled:
%{
figure(2)
plot(stresstolookat,z)
hold on
plot(culledstresses,z(z>(max(z)-0.5*(max(z)-min(z)))))
hold off
%}

%% MULTISTART 
try
    % Now seed based on the culled stresses, 4 points to find the minima. This
    % is fairly efficient, and mostly robust. 
    seedpoints=unique(round(linspace(1,numel(culledstresses),4)));

    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective',...
        fun,'x0',minmu+0.5*rangmu,'lb',minmu,'ub',maxmu,'options',opts);
    ms = MultiStart;
    ms.Display='off';
    tpoints = CustomStartPointSet(culledstresses(seedpoints)');

    [beststress,f] = run(ms,problem,tpoints);
catch
    %If you do not have createOptimproblem
    stressestolookat=minmu:5:maxmu; %Finer mesh

    z=stressestolookat;
    for i=1:numel(stressestolookat)
        z(i)=-fun(stressestolookat(i));
    end
    if nnz(isnan(z))>0
        error('Nans in the utility function')
    end
    I=find(z==max(z(:)));
    if nnz(I)>1
        beststress=min(stressestolookat(I));
    else
        beststress=stressestolookat(I(1));
    end
    f=z(I(1));
end

%% Debug scripts
%LOOK AT THE UTILITY WITHIN A RANGE
%{
    stresstolookat=100:5:700;
    z=stresstolookat;
    for i=1:numel(stresstolookat)
        z(i)=-fun(stresstolookat(i));
    end
    % figure();
%     z=z-min(z);%NORMALISE Z FOR CONVENIENCE
%     z=z./max(z);
    f=max(z);
    figure(3)
    plot(stresstolookat,z)
    xlabel('Test Stress (MPa)');
    ylabel('Utility Function')
    hold on
    scatter(beststress,f)
    pause(0.01)
%}

% PLOT THE PAPER FIGURE WITH THE VARIATION IN UTILITY FOR SOME NUMBER OF
% SPECIMENS
%{
%Figure out which sample you're in
k=evalin('caller','k');%CMM HACK : copy this line into B_simulate as well in order for this to work. 

whichsamp=[1,2,3,20,60,90]; %which samples do you want to look at 
if ismember(k,whichsamp)
    stresstolookat=100:5:700;
    z=stresstolookat;
    for i=1:numel(stresstolookat)
        z(i)=-fun(stresstolookat(i));
    end
    figure(1);
    z=z-min(z);%NORMALISE Z FOR CONVENIENCE
    z=z./max(z);
    z=z-(1+find(k==whichsamp));
    plot(stresstolookat,z);
    xlabel('Test Stress');
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
        print(['utilityfirst',num2str(whichsamp(end))],'-dpng','-r0')
        savefig(['utilityfirst',num2str(whichsamp(end)),'.fig'])
    end
end
%}
end