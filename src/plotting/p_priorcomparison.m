function p_priorcomparison(lpriorstack,theta,sigma,leg,xl,yl)
%compare between priors showing marginal mean and standard deviations as
%well as joint posterior

%INPUTS: cell array of priors to be compared, theta, sigma, legend name cell
%array, sigma limits for plot, theta limits for plot

subplot(2,2,3)
co='grbcmyrgbcmyk'; %colors for plotting
for i=1:size(lpriorstack,1)
    [C,h(i)]=p_contourHPD(lpriorstack{i},'newfig',0); %create the HPD contours
    %fill in the contours
    h(i).LineColor=co(i);
    c2=C;
    I=find(c2(1,:)<1);
    c2=c2(:,I(2)+1:end);
    x1=c2(1,:);
    y1=c2(2,:);
    C(:,C(1,:)<1)=[];
    x2 = C(1,:);
    y2 = C(2,:);
    hold on
    fill(x2,y2,co(i),'FaceAlpha',0.3,'LineStyle','none')
    fill(x1,y1,co(i),'FaceAlpha',0.6,'LineStyle','none')
end
ylim(yl)
xlim(xl)
legend(h(1:i),leg); %use string to fill in legend
ylabel('\theta /MPa')
xlabel('\sigma /MPa')
%create marginal probability plots off to the side
subplot(8,2,7)
linew=1.5;
for i=1:size(lpriorstack,1)
    plot(sigma,sum(exp(lpriorstack{i}),1),co(i),'LineWidth',linew)
    hold on
end
set(gca,'xticklabel',[])
ylabel('P_\sigma')
xlim(xl)
%create marginal probability plots off to the side
subplot(2,8,13)
for i=1:size(lpriorstack,1)
    plot(sum(exp(lpriorstack{i}),2),theta,co(i),'LineWidth',linew)
    hold on
end
ylim(yl)
set(gca,'yticklabel',[])
xlabel('P_\mu')
end