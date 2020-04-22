function [ vf ] = binning( vx,vpts, wind )
% 
% Description
% 
% Performs a binning (histogram estimator) of the data on vx, and
% evaluates the estimated density function on the points on vpts. The
% length of the window used is fixed on wind.
% 
% Inputs
% vx - response data that will be used to perform the binning
% vpts - explanatory variable used to perform the binning of vx
% wind - length of the windom used in the binning
% 
% Outputs
% vf - binned response varible

n = length(vx);
N = length(vpts);
vf = zeros(N,1);

for ii=1:N
    fhat = sum((vpts(ii)-(0.5*wind)<=vx).*(vx<=vpts(ii)+(0.5*wind)))/(wind*n);
    vf(ii) = fhat;
end

end

