function y = IDWT(cA, cD, filter)
% IDWT
% Single-level Inverse Discrete Wavelet Transform.
%
% Syntax
%   y = IDWT(cA, cD, wavfilter) 
%
% Description
% IDWT calculates the single-level reconstructed data y with respect to the
%   smoothed coefficients cA, the detail coefficients cD and the wavelet
%   filter 'filter'.
%
% Inputs
% cA            - Smoothed (scaling) coefficients.
% cD            - Detail coefficients.
% filter        - Scaling wavelet filter.
%
% Outputs
% y             - Reconstructed data.
%
% Example
% x = [1 0 -3 2 1 0 1 2];
% wfilter = [sqrt(2)/2 sqrt(2)/2];
% [C, D] = DWT(x, wfilter);
% y = IDWT(C, D, wfilter)
%
% % Compare the reconstructed data y with the original
% % data x by the sum of the absolute values of the
% % diferences
% sum(abs(x-y))

% By Michel H. Montoril at Georgia Institute of Technology, Sep 2014.


    N = length(filter);
    n = 2*length(cA);
    
    Lo_R = filter;
    Hi_R = filter(N:-1:1);
    Hi_R(2:2:N) = -Hi_R(2:2:N);
    
    w  = mod(0:(N/2-1), n/2) + 1;              % Make periodic
    
    C(1:2:(n+N)) = [cA cA(1,w)];               % Upsample & keep periodic
    D(1:2:(n+N)) = [cD cD(1,w)];               % Upsample & keep periodic
    y = conv(C, Lo_R) + conv(D, Hi_R);         % Convolve & add
    y = y((N:N+n-1)-1);                        % Periodic part
    
end
