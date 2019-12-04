function y = wr(wc, J0, filter)
% wr
% Wavelet reconstruction of a signal by Inverse Discrete Wavelet Transform.
%
% Syntax
%   y = wd(wc, J0, filter) 
%
% Description
% wr reconstructs a signal that is in wavelet domain
% with a specific coearsest level J0.
%
% Inputs
% wc            - Wavelet coefficients to be reconstructed.
% J0            - Coarsest level of the decomposition.
% filter        - Scaling wavelet filter.
%
% Outputs
% y             - The signal reconstructed from the wavelet coefficients.
%
%
% Example
% % Let us reconstruc the signal x below
% x = [1 0 -3 2 1 0 1 2];
%
% % First we decompose it using the wd function with
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
%
% % The reconstructed signal is obtaind by
% y = wr(xwd, j0, wfilter);
% % The sum of the absolute deviations between x and y is
% sum(abs(x-y))

% By Michel H. Montoril at Georgia Institute of Technology, Sep 2014.

    n = length(wc);
    J = log2(n);
        
    scales = J0:J;
    sqrscales = 2.^scales;
    
    if J == J0
        y = wc;
        return
    end
    
    y = IDWT(wc(1:sqrscales(1)), wc((sqrscales(1)+1):sqrscales(2)), filter);
    
    for i = 2:(J-J0)
        y = IDWT(y, wc((sqrscales(i)+1):sqrscales(i+1)), filter);
    end

end
