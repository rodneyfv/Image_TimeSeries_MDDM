function [ fEfr1 ] = mixtureProbs(y, x, s, delt, min_pts, wJ, wfilt, wprec)
%
% Description
% This function is used to apply the method of Montoril et al (2019) in the
% vector y.
%
% Inputs
% y             - Response variable (not in the regression context)
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
%
% Outputs
% fEfr1       - function with the probability of each point belong to the
%                first group

% threshold used to separate the observations in y in two groups
tmp = quantile(y,0.9);
% estimates of the mean in the two groups separated above
muG = mean(y(y>tmp));
muL = mean(y(y<=tmp));

% Applying the transformation to take the problem to a regression context
w = (y - muL)/(muG - muL);

% Obtaining the raw scaling coefficient estimates
% We need to use a row vector in the function below
[w1, ~] = WavCoefEsts(w', x, s, delt, min_pts, wJ, wfilt, wprec);

% Decomposing the raw estimates into wavelet domain
[wd1, snoise1] = wd(w1, 0, wfilt, true);

% Shrinkage using Efromovich's approach
[rwdEfr1, Jest1, ~] = WavShrink(wd1, wfilt, true, 'Efromovich', true, snoise1);

% Function estimates regularized by Efromovich's shrinkage
yy1 = PHI(x, Jest1, wfilt, wprec, true);
fEfr1 = yy1*rwdEfr1';
fEfr1(fEfr1<0) = 0;fEfr1(fEfr1>1) = 1;

end

