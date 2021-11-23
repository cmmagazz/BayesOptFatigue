
figure
plot(ResultSet.raw.failurestress(:,1)-(ResultSet.raw.failurestress(:,2)-3)*ResultSet.details.step.stepsize,'k-')
xlabel('Sample Number')
ylabel('Starting Stress /MPa')
% 
% print(['startingstress'],'-dpng','-r0')
% savefig(['startingstress.fig'])
%%

% whichsamp=[2,3,15,25,45,70]; %which samples do you want to look at 
whichsamp=[70]; %which samples do you want to look at 

minstressstep=ResultSet.details.step.stepsize;
startingstress=ResultSet.details.startingstress;

startingstresses=linspace(startingstress-minstressstep,startingstress,50);

numstresssteptolookat=1;
if numstresssteptolookat ==1
    stepsizetolookat=minstressstep;
else
    stepsizetolookat=linspace(minstressstep,100,numstresssteptolookat);
end

for j=whichsamp
    ResultSet2=ResultSet;
    if j==1
        ResultSet2.raw.failurestress=[];
        ResultSet2.raw.lprior=g_calcprior([],ResultSet.details.theta);

    else
        ResultSet2.raw.failurestress=ResultSet2.raw.failurestress(1:j-1,:);
        ResultSet2.raw.lprior=g_calcprior(ResultSet2);

    end
    fun=@(x)g_STEP_UTILITY(x,ResultSet2);

    z=zeros(1,numel(startingstresses));
    second=zeros(1,numel(startingstresses));
    for i=1:numel(startingstresses)
        [z(i),second(i)]=fun(startingstresses(i));
    end
    z=-1.*z;
    figure(9);
    z=z-min(z);%NORMALISE Z FOR CONVENIENCE
    z=z./max(z);
    z=z-(1+find(j==whichsamp));
    plot(startingstresses+320,z);
    xlabel('Starting Stress');
    ylabel('Utility Function')
    hold on
        I=find(z==max(z(:)));
    beststress=min(startingstresses(I));

    scatter(beststress+320,max(z),'kx')
    if j==whichsamp(end)
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

%%

numstresssteptolookat=50;

startingstresses=linspace(startingstress-minstressstep,startingstress,40);

stepsizetolookat=linspace(10,160,numstresssteptolookat);

startingstresses=0:0.2:10;%<<<<<<<<

NS1=combvec(stepsizetolookat,startingstresses);

z=zeros(1,size(NS1,2));
second=zeros(1,size(NS1,2));
w = waitbar(0);


for i=1:numel(NS1(1,:))
    waitbar(i/numel(NS1(1,:)),w,['Sample number: ' num2str(i),'/',num2str(numel(NS1(1,:)))])
    ResultSet.details.step.stepsize=NS1(1,i);
    ResultSet.details.startingstress=startingstress-NS1(2,i).*NS1(1,i)./10;%<<<<<<<<<<<<<<
    fun=@(x)g_STEP_UTILITY(x,ResultSet);
%     [z(i),second(i)]=fun(NS1(2,i));
    [z(i),second(i)]=fun(startingstress-NS1(2,i).*NS1(1,i)./10);
end
close(w)
z=-1.*z;
if nnz(isnan(z))>0
    error('Nans in the utility function')
end


    gridN=reshape(NS1(1,:),numel(stepsizetolookat), numel(startingstresses));
    gridS=reshape(NS1(2,:),numel(stepsizetolookat), numel(startingstresses));
    gridz=reshape(z,numel(stepsizetolookat), numel(startingstresses));
    h=contourf(gridN,gridS*10,gridz,1000,'LineColor','None');
    xlabel('Stress Step /MPa');
    ylabel('Starting Stress Pertubation /%')
    c=colorbar;
    c.Label.String = 'Utility';
    
    %%
print(['utility2d'],'-dpng','-r900')
        savefig(['utility2d.fig'])
        
        
        save('steputilitydata.mat')