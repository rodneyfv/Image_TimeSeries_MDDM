function [ y ] = normalize( x )
%
% Description
% Function used to normalize a vector
%
% Inputs
% x          - vector nx1 with data
% Outputs
% y          - vector nx1 with mean one and variance 1

y = (x - mean(x))./std(x);

end

