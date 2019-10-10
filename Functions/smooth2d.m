function [ mHs ] = smooth2d( mH, lambda )
% Function that applies to a matrix H the smoothing method 
% described in Eilers and Goeman (2004, Bioinformatics).
% First, columns are smoothed, then rows are smoothed.
% lambda is a tunning parameter.

% 2-D exponential smoother
if length(lambda) == 1
    lambda = [lambda, lambda];
end

% Here, index x is related to variables used in the column's
% smoothing, and index y are variables related to the row's
% smoothing

% Build penalty matrices
[m,n] = size(mH);
% identity matrices
mEx = eye(m);
mEy = eye(n);
% matrices of first differences
mD1x = diff(mEx);
mD1y = diff(mEy);
% matrices of second differences
mD2x = diff(mD1x);
mD2y = diff(mD1y);
% matrices described in Eq.(8) of Eilers and Goeman (2004)
mQx = mEx + (lambda(1)^2)*(mD2x' * mD2x) + (2*lambda(1)) * (mD1x'*mD1x);
mQy = mEy + (lambda(1)^2)*(mD2y' * mD2y) + (2*lambda(1)) * (mD1y'*mD1y);
% First we smooth the columns by solving the system 
% Qx*Z.hat=H for Z.hat, then we take the transpose of Z.hat,
% and smooth the rows, which correspond to solve the system
% Qy*Z2.hat=t(Z.hat). Our smoothed H then is given by Z2.hat
mHs = linsolve(mQy,linsolve(mQx,mH)')';

end

