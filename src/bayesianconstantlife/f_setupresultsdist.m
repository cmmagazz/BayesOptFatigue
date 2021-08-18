function ResultSetdetails=f_setupresultsdist(ranges,dist,NUMEL)
%Function to setup the cell of vectors of model parameter values
%   For Normal:
%       dist='norm';
%       mu = mean fatigue strength
%       sigma = standard deviation
%   For 2-param Weibull: 

if isempty(dist) %set auto value
    dist = 'norm';
end
if strcmp(dist,'norm') || strcmp(dist,'lognorm') || strcmp(dist,'2pwbl') || strcmp(dist,'type1')
    needsize=[2,2];
    if nnz(size(ranges)~=needsize)>0
        error('Need 2x2 for ranges')
    end
    minmu=ranges(1,1);
    maxmu=ranges(1,2);
    mu=linspace(minmu,maxmu,NUMEL);%from min to max in steps of

    minsigma=ranges(2,1);
    maxsigma=ranges(2,2);
    sigma=linspace(minsigma,maxsigma,NUMEL);

    %Insert into ResultSet struct
    ResultSetdetails.theta={mu, sigma};
    ResultSetdetails.numparams=2;
elseif strcmp(dist,'3pwbl') || strcmp(dist,'gev')
    needsize=[3,2];
    if nnz(size(ranges)~=needsize)>0
        error('Need 3x2 for ranges')
    end
    minA=ranges(1,1);
    maxA=ranges(1,2);
    A=linspace(minA,maxA,NUMEL); %

    minB=ranges(2,1);
    maxB=ranges(2,2);
    B=linspace(minB,maxB,NUMEL);    
    
    minC=ranges(3,1);
    maxC=ranges(3,2);
    C=linspace(minC,maxC,NUMEL);

    %Insert into ResultSet struct
    ResultSetdetails.theta={A, B, C};
    ResultSetdetails.numparams=3;
end
ResultSetdetails.dist=dist;
end