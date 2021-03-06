function prob=g_calcprobCDF(stress,theta,dist)
% More generalised form to calculate the probability (CDF) of recording some result. Supported distributions: 
%   - Normal                = 'norm'
%   - Log Normal            = 'lognorm'
%   - 2 Parameter Weibull   = '2pwbl'
%   - 3 Parameter Weibull   = '3pwbl'
%   - Generalised Extreme Value: 
%       - Frechet           = 'gev'
%       - Gumbel            = 'type1'

if strcmp(dist,'norm')
    [A,B]=meshgrid(theta{1},theta{2});
    prob=normcdf(stress,A',B');
elseif strcmp(dist,'lognorm')
    [A,B]=meshgrid(theta{1},theta{2});
    prob=logncdf(stress,log(A'),B');
elseif strcmp(dist,'2pwbl')
    [A,B]=meshgrid(theta{1},theta{2});
    prob=wblcdf(stress,A',B');
elseif strcmp(dist,'3pwbl')
    [A,B,C]=meshgrid(theta{1},theta{2},theta{3});
    A=permute(A,[2,1,3]);
    B=permute(B,[2,1,3]);
    C=permute(C,[2,1,3]);
%     threepwblcdf = @(x,a,b,c) (x>c).*(1-exp(-((x-c)./a).^b));
%     prob=threepwblcdf(stress,A,B,C);
    wherelarger=(stress-C)>0;%deal with requirement for non-imaginary outputs
    
    if ndims(A)>2
        prob=zeros(size(A));
        prob(wherelarger) = (1-exp(-((stress-C(wherelarger))./A(wherelarger)).^B(wherelarger)));
    elseif isvector(stress) || ismatrix(stress)
        prob=zeros(size(stress));
        prob(wherelarger) = (1-exp(-((stress(wherelarger)-C)./A).^B));
    end
elseif strcmp(dist,'gev')
    [A,B,C]=meshgrid(theta{1},theta{2},theta{3});
    A=permute(A,[2,1,3]);
    B=permute(B,[2,1,3]);
    C=permute(C,[2,1,3]);
    prob = gevcdf(stress,C,B,A);
elseif strcmp(dist,'type1')
    [A,B]=meshgrid(theta{1},theta{2});
    prob = gevcdf(stress,zeros(size(A')),B',A');
end

end
