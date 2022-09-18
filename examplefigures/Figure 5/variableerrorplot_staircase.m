%% A tool for the experimentalist to view error variability across multiple number of samples
%Figure 5c) and d) cuts the curve of Figure 5a) and b) at a certain number
%of samples, to see how - as a function of step size - the error changes.
%The number of samples chosen in the paper is 100, representing the end of
%the simulation. If you want to view the cut at any number of samples, or
%at a range, please use the code below, changing the variable
%"howmanysamples"


for howmanysamples=5:5:100
    errorvariationplot(howmanysamples,progthet, progsig, distwidth, protocolstotest, widths)
    drawnow()
    clf
end

function errorvariationplot(howmanysamp, progthet, progsig, distwidth, protocolstotest, widths)
errbarthet=(squeeze(progthet(:,:,:,:,1))-400).*(100/400);
errbarsig=(squeeze(progsig(:,:,:,:,1))-distwidth).*(100/distwidth);

colours={'r','g'};

%First look at mean
arraytotest=sqrt(errbarthet(:,:,:,:).^2);
%determine error value and stdev
for protocolk=1:length(protocolstotest)
    for widthl=1:numel(widths)
        erroratN(protocolk,widthl)=mean(arraytotest(protocolk,widthl,:,howmanysamp));
        errbar(protocolk,widthl)=std(arraytotest(protocolk,widthl,:,howmanysamp));
    end
end
%plot using shadederrorbar
for protocolk=1:length(protocolstotest)
    subplot(1,2,1)
    shadedErrorBar(widths,erroratN(protocolk,:),errbar(protocolk,:),'lineProps',colours{protocolk});
    h=gca;
    set(h,'xscale','log')
end


xlabel('Ratio of step size to \sigma')
ylabel('Error in Mean, %')
ylim([0 inf])

% AGAIN FOR SIGMA
arraytotest=sqrt(errbarsig(:,:,:,:).^2);
for protocolk=1:length(protocolstotest)
    for widthl=1:numel(widths)
        erroratN(protocolk,widthl)=mean(arraytotest(protocolk,widthl,:,howmanysamp));
        errbar(protocolk,widthl)=std(arraytotest(protocolk,widthl,:,howmanysamp));
    end
end
for protocolk=1:length(protocolstotest)
    subplot(1,2,2)
    shadedErrorBar(widths,erroratN(protocolk,:),errbar(protocolk,:),'lineProps',colours{protocolk});
    h=gca;
    set(h,'xscale','log')
end

legend({'Staircase','Bayesian Staircase'})
xlabel('Ratio of step size to \sigma')
ylabel('Error in Standard Deviation, %')
ylim([0 inf])

sgtitle(['Convergence after ' num2str(howmanysamp) ' samples'],'FontSize',12)

end
