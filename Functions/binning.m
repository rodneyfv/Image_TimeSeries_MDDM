function [ vf ] = binning( vx,vpts, wind )
% Performs a binning (histogram estimator) of the data on vx, and
% evaluates the estimated density function on the points on vpts. The
% length of the window used is fixed on wind.

n = length(vx);
N = length(vpts);
vf = zeros(N,1);

for ii=1:N
    fhat = sum((vpts(ii)-(0.5*wind)<=vx).*(vx<=vpts(ii)+(0.5*wind)))/(wind*n);
    vf(ii) = fhat;
end

end

