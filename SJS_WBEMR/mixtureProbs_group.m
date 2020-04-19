function [ fHat ] = mixtureProbs_group(y, yg, x, s, delt, min_pts, wJ, wfilt, wprec, rawest, estimator)
%
% Description
% This function is used to apply the method of Montoril et al (2019) in the
% vector y. The k-means algorithm is applied to choose two different groups
% whose means are computed and used in Montoril's method.
%
% Inputs
% y             - Response variable (not in the regression context)
% yg            - Pre-separation of variables in two groups (0s and 1s)
% x             - Regressor variable in the regression context.
% s             - Value used to define the intervals $ [x_{i-s},x_{i+s}] $.
% delt      - Length of the sub-intervals to be considered in the
%                   Trapezoidal rule.
% min_points          - Auxiliary argument. The number of points to be used in
%                   the Trapezoidal rule, in case intDelta is not small
%                   enough.
% wJ          - Resolution level of approximation.
% wfilt     - Scaling wavelet filter.
% wprec       - Precision approximation measured by the number of
%                   Daubechies-Lagarias steps.
% rawest      - Type of raw estimator. If 'wavcoef', Eq. (8) of Montoril et
%           al (2019) is used, or if 'wavcoefint', their Eq. (9) is used.
% estimator   - Type off estimator, which can be 'HardThresh',
%           'Efromovich', or 'LocLinReg'.
% Outputs
% fHat       - function with the probability of each point belong to the
%                first group

% estimates of the mean in the two groups separated above
muG = mean(y(yg==1));
muL = mean(y(yg==0));
if(muG < muL) % making sure muG is the largest mean
    tmp = muG;
    muG = muL;
    muL = tmp;
end

% Applying the transformation to take the problem to a regression context
w = (y - muL)/(muG - muL);

% Obtaining the raw scaling coefficient estimates. If the user sets rawest
% as 'wavcoefint', the raw estimator of Montoril et al (2019) that uses
% numerical integration is employed, otherwise the method presented in 
% their equation (8) is used.
if(strcmp(rawest,'wavcoefint'))
    % We need to use a row vector in the function below.
    % Raw estimator of equation (9) in Montoril et al (2019)
    [~, w1] = WavCoefEsts(w', x, s, delt, min_pts, wJ, wfilt, wprec);
elseif(strcmp(rawest,'wavcoef'))
    % Raw estimator of equation (8) in Montoril et al (2019)
    [w1,~] = WavCoefEsts(w', x, s, delt, min_pts, wJ, wfilt, wprec);
else
    error('Raw estimator not defined correctly')
end

% Decomposing the raw estimates into wavelet domain
[wd1, snoise1] = wd(w1, 0, wfilt, true);

if strcmp(estimator,'HardThresh')
    % Decomposing the raw estimates into wavelet domain
    [wd1, snoise1] = wd(w1, 0, wfilt, true);
    wdH1 = wd1(1:2^wJ);
    % Obtaining the raw function estimates
    yy = PHI(x, wJ, wfilt, wprec, true);
    % Shrinkage using hard threshold
    rwdH1 = WavShrink(wdH1, wfilt, true, 'Hard', true, snoise1);
    % Function estimates regularized by Hard threshold
    fHat = yy*rwdH1';
    fHat(fHat<0) = 0;
    fHat(fHat>1) = 1;
    
elseif strcmp(estimator,'Efromovich')
    % Shrinkage using Efromovich's approach
    [rwdEfr1, Jest1, ~] = WavShrink(wd1, wfilt, true, 'Efromovich', true, snoise1);
    % Function estimates regularized by Efromovich's shrinkage
    yy1 = PHI(x, Jest1, wfilt, wprec, true);
    fHat = yy1*rwdEfr1';
    fHat(fHat<0) = 0; fHat(fHat>1) = 1;  
    
elseif strcmp(estimator,'LocLinReg')
    % Obtaining the raw function estimates
    yy = PHI(x, wJ, wfilt, wprec, true);
    f1row1 = yy*w1';
    % Function estimates regularized by Local Linear Regression
    % bandwidth used in the ksrlin function
    r.n=length(x);
    hx=median(abs(x-median(x)))/0.6745*(4/3/r.n)^0.2;
    hy=median(abs(f1row1-median(f1row1)))/0.6745*(4/3/r.n)^0.2;
    h=sqrt(hy*hx);
    % the fourth argument in ksrlin indicates we want the results in n
    % points
    r1=ksrlin(x,f1row1,h,r.n);
    r1.f(r1.f<0) = 0;
    r1.f(r1.f>1) = 1;
    fHat = r1.f;
    
else
    error('Estimation method not defined correctly')
end
    
end

