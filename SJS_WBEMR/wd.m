function [y, sig] = wd(signal, J0, filter, usemad)
% wd
% Wavelet decoposition of a signal by Discrete Wavelet Transform.
%
% Syntax
%   [y, sig] = wd(signal, J0, filter, usemad) 
%
% Description
% wd decomposes a periodic signal into wavelet domain from
% the finest level to the a specified coarsest level J0.
%
% Inputs
% signal        - Data set to be decomposed.
% J0            - Coarsest level in the decomposition.
% filter        - Scaling wavelet filter.
% usemad        - Logical argument. If it is true, the function will
%                   calculate the median absolute deviation (MAD) of
%                   the detail coefficients from the finest level.
%                   Otherwise, the function will calculate the sample
%                   standard deviation of the detail coefficients from
%                   the finest level.
%
% Outputs
% y             - The signal decomposed into wavelet domain.
% sig           - If the signal has noise, this output will provide an
%                   estimate of the standard deviation of the noise,
%                   based on the finest detail coefficients. This
%                   estimate may be based on the MAD or the standard
%                   deviation. See the usemed input for details.
%
% Example
% % Let us consider the signal x below
% x = [1 0 -3 2 1 0 1 2];
% % the a 10-tap Daubechies Least Asymetric filter called wfilter
% wfilter = [0.027333068345164, 0.029519490926073,-0.039134249302583,...
%            0.199397533976996, 0.723407690403808, 0.633978963456949,...
%            0.016602105764424,-0.175328089908107,-0.021101834024930,...
%            0.019538882735387];
% % and the coarsest resolution level j0
% j0 = 0;
% % The wavelet decomposition of x (xwd) and the MAD of the finest
% % detail coefficients (dsig) are obtained by
% [xwd, dsig] = wd(x, j0, wfilter, true);

% By Michel H. Montoril at Georgia Institute of Technology, Sep 2014.

    n = length(signal);
    J = log2(n);
    scales = (J-1):-1:J0;
    ncols = 2.^scales;
    cncols = cumsum(ncols);
    y = zeros(1, n);
    [C, D] = DWT(signal, filter);
    
    if usemad==true
        sig = 1.4826*mad(D, 1);
    else
        sig = std(D);
    end
    
    y(1:ncols(1)) = D(end:-1:1);

    for i = 2:(J-J0)
        [C, D] = DWT(C, filter);
        y((cncols(i-1)+1):cncols(i)) = D(end:-1:1);
    end
    
    y = y(end:-1:1);
    y(1:2^J0) = C;
    
end