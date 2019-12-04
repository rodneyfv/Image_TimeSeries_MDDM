function [ca, cd] = DWT(x, filter)
% DWT
% Single-level Discrete Wavelet Transform.
%
% Syntax
%   [ca, cd] = DWT(x, wavfilter) 
%
% Description
% DWT calculates the smoothed coefficients ca and the detail coefficients
%   cd of a periodic data x for a specific scaling wavelet filter 'filter'.
% 
% Inputs
% x             - Data set to be decomposed.
% filter        - Scaling wavelet filter.
%
% Outputs
% ca            - Smoothed (scaling) coefficients.
% cd            - Detail coefficients.
%
% Example
% x = [1 0 -3 2 1 0 1 2];
% wfilter = [sqrt(2)/2 sqrt(2)/2];
% [C, D] = DWT(x, wfilter)

% By Michel H. Montoril at Georgia Institute of Technology, Sep 2014.

    N = length(filter);
    n = length(x);
    
    Lo_D = filter(end:-1:1);
    Hi_D = filter;
    Hi_D(1:2:end) = -Hi_D(1:2:end);
    
    ca = [x(mod(((1 - N):-1), n) + 1)  x]; % make periodic
    cd = conv(ca, Hi_D);                   % Convolve,
    cd = cd((N:2:(N + n - 2)) + 1);          % keep periodic and decimate
    ca = conv(ca, Lo_D);                   % Convolve,
    ca = ca((N:2:(N + n - 2)) + 1);          % keep periodic and decimate
    
end
