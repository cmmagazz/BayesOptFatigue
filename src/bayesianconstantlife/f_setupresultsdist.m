function ResultSetdetails=f_setupresultsdist(ranges,dist,NUMEL)
%Function to setup the cell of vectors of model parameter values
%   For Normal:
%       dist='norm';
%       mu = mean fatigue strength
%       sigma = standard deviation
%   For 2-param Weibull: 


if strcmp(dist,'norm') || isempty(dist)
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

    ResultSetdetails.dist='norm';
elseif strcmp(dist,'lognorm') 
    needsize=[2,2];
    if nnz(size(ranges)~=needsize)>0
        error('Need 2x2 for ranges')
    end
    minmu=ranges(1,1);
    maxmu=ranges(1,2);
    mu=linspace(minmu,maxmu,NUMEL);
%     mu=minmu:2:maxmu; %from min to max in steps of 2

    minsigma=ranges(2,1);
    maxsigma=ranges(2,2);
    sigma=linspace(minsigma,maxsigma,NUMEL);

    %Insert into ResultSet struct
    ResultSetdetails.theta={mu, sigma};

    ResultSetdetails.dist='lognorm';
elseif strcmp(dist,'2pwbl') 
    needsize=[2,2];
    if nnz(size(ranges)~=needsize)>0
        error('Need 2x2 for ranges')
    end
    minA=ranges(1,1);
    maxA=ranges(1,2);
    A=linspace(minA,maxA,NUMEL); %from min to max in steps of 2

    minB=ranges(2,1);
    maxB=ranges(2,2);
    B=linspace(minB,maxB,NUMEL);

    %Insert into ResultSet struct
    ResultSetdetails.theta={A, B};

    ResultSetdetails.dist='2pwbl';
elseif strcmp(dist,'3pwbl') 
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

    ResultSetdetails.dist='3pwbl';
end
end