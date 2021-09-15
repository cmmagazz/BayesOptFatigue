toplot=[1,10,30,60];
ResultSet.details=f_setupresultsdist([[250, 550];[1, 150]],'norm',300);

for i=1:numel(toplot)
    [lprior,~,shannon]=g_calcprior(ResultSet.raw.failurestress(1:toplot(i),:),ResultSet.details.theta);
    subplot(2,2,i)
    p_HPD(lprior,'newfig',0)
    xticks(0:20:150)
    if i==1
        title([num2str(toplot(i)),' Sample'])
    else 
        title([num2str(toplot(i)),' Samples'])
    end
    pause(1)
    
end


print('p_HPD evolution','-dpng','-r0')
savefig(['p_HPD evolution','.fig'])