function prob=g_calcprobCDF(stress,theta,dist)
if strcmp(dist,'norm')
    prob=normcdf(stress,repmat(theta{1}',[1,size(theta{2},2)]),repmat(theta{2},[size(theta{1},2),1]));
elseif strcmp(dist,'lognorm')
    prob=logncdf(stress,repmat(theta{1}',[1,size(theta{2},2)]),repmat(theta{2},[size(theta{1},2),1]));
elseif strcmp(dist,'2pwbl')
    prob=wblcdf(stress,repmat(theta{1}',[1,size(theta{2},2)]),repmat(theta{2},[size(theta{1},2),1]));
elseif strcmp(dist,'3pwbl')
    [A,B,C]=meshgrid(theta{1},theta{2},theta{3});
    threepwblcdf = @(x,a,b,c) (x>c).*(1-exp(-((x-c)./a).^b));
    prob=threepwblcdf(stress,A,B,C);
end

end
