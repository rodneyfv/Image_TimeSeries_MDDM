function [wavcoef, wavcoefint] = WavCoefEsts(y, x, sn, intDelta, intN, intJ, intFilter, intPrec)
% WavCoefEsts
% Wavelet coefficients estimated in the context of heteroscedastic
%   nonparametric regression model.
%
% Syntax
%   [wavcoef, wavcoefint] = WavCoefEsts(y, x, sn, intDelta, intN, intJ, intFilter, intPrec) 
%
% Description
% Wavelet coefficients estimated in the context of the raw estimators in
%   Section 3 of the paper. Basically, the command returns two different
%   wavelet coefficient estimates of the function $ f $ in a
%   heteroscedastic regression model y = f(x) + e(f,x), where $ e(x, f)
%   corresponds to the error and it depends on x and the function f.
%
% Inputs
% y             - Response variable in the regression context.
% x             - Regressor variable in the regression context.
% s             - Value used to define the intervals $ [x_{i-s},x_{i+s}] $.
% intDelta      - Length of the sub-intervals to be considered in the
%                   Trapezoidal rule.
% intN          - Auxiliary argument. The number of points to be used in
%                   the Trapezoidal rule, in case intDelta is not small
%                   enough.
% intJ          - Resolution level of approximation.
% intFilter     - Scaling wavelet filter.
% intPrec       - Precision approximation measured by the number of
%                   Daubechies-Lagarias steps.
%
% Outputs
% wavcoef       - Wavelet coefficient estimates based on the first raw
%                   estimator in the paper.
% wavcoefint    - Wavelet coefficient estimates based on the second raw
%                   estimator in the paper.
%
% Details
% Observe that the argument s used here is slightly different from the used
% in our paper. In the paper, s define the amplitue of the intervals
% $ [x_{i-s/2},x_{i+s/2}] $, and here the intervals are
% $ [x_{i-s},x_{i+s}] $. Thus, s=1 here is equivalent to s=2 in our paper.
%
% Example
% % See its use in the Application.m file

% By Michel H. Montoril at Georgia Institute of Technology, Sep 2014.

    %#codegen
    coder.inline('never')
    
    [phi_mat, ~] = PHI(x, intJ, intFilter, intPrec, true);
    
    lowbound = [zeros(1, sn) x(1:(end - sn))];
    upbound = [x((sn + 1):end) ones(1, sn)];
    
    wavcoef = (1/(2*sn)) * (upbound-lowbound) .* y * phi_mat;
    
    intMatrix = WavIntegral(lowbound, upbound, intDelta, intN, intJ, intFilter, intPrec);%, phi_k(1));
    
    wavcoefint = (1/(2*sn)) * y * intMatrix;
    
    
end
%------------------------------------------------------%
%-  This function calculetes the integral of wavelets -%
%------------------------------------------------------%
function wavint = WavIntegral(lowerbound, upperbound, delta, nint, J, filter, prec)
    n = length(lowerbound);
    
    p = 2^J;
    N = length(filter);
    kmin = 1 - N;
    kmax = p - 1;
    wavint = zeros(n, kmax - kmin + 1);
    
    for i = 1:n
        leng = max( nint, ceil( (upperbound(i) - lowerbound(i))/delta ) );
        xint = linspace(lowerbound(i), upperbound(i), leng);
        [Yint, kvals] = PHI(xint, J, filter, prec, false);
        wavint(i,(kvals(1)-kmin+1):(kvals(2)-kmin+1)) = Trapezoidal(Yint, xint(2)-xint(1));
    end
    
    wavint = sum(reshape([wavint zeros(n, p - mod(size(wavint, 2), p))], n, p, []), 3);
    wavint = wavint(:, mod(1:p, p) + 1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  Functions needed   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------%
%-  This function calculetes the numerical integral   -%
%-  of matrices by the trapezoidal rule               -%
%------------------------------------------------------%
function res = Trapezoidal(Y, delta)
    res = delta * sum(Y) - (delta/2) * (Y(1,:) + Y(end,:));
end

% %------------------------------------------------------%
% %-  This function calculetes the design matrix        -%
% %-  of a specific wavelet basis by the                -%
% %-  Daubechies-Lagarias alorithm                      -%
% %------------------------------------------------------%
% function w = MPHI(x, j, filt, kmin, kmax, prec, periodized)
%    n = length(x);
%    w = zeros(n, kmax - kmin + 1);
%    for i = 1:n
%        w(i,:) = Phi(x(i), j, filt, prec, kmin, kmax);
%    end
%    
%    if periodized == true
%        w = sum(reshape([w zeros(n, 2^j - mod(size(w, 2), 2^j))], n, 2^j, []), 3);
%    end
% end
% 
% %------------------------------------------------------%
% %- This function calculates the wavelet basis in a    -%
% %- specific point by de Daubechies-Lagarias algorithm -%
% %------------------------------------------------------%
% function phi = Phi(x, filter, N, prec, kmin, kmax)
%     int = floor(x);
%     dec = x - int;
%     dy = dec2bin(dec, prec);
%     k1 = ceil(x - N + 1);
%     k2 = floor(x - 1e-9);
%     T0 = Tmat(filter, N, 0);
%     T1 = Tmat(filter, N, 1);
%     prod = eye(N - 1);
%     for i = 1:prec
%         if dy(i)==1
%             prod = prod*T1;
%         else
%             prod = prod*T0;
%         end
%     end
%     phi = zeros(1, kmax - kmin + 1);
%     phi((k2-kmin+1):-1:(k1-kmin+1)) = mean(prod, 2);
% end
% 
% %------------------------------------------------------%
% %- This function calculates the first n elements of a -%
% %- diadic representation of x                         -%
% %------------------------------------------------------%
% function a = dec2bin(x, n)
%     a = zeros(1, n);
%     y = x;
%     for i = 1:n
%          y = 2*y;
%         a(i) = floor(y);
%         y = y - floor(y);
%     end
% end
% 
% %------------------------------------------------------%
% %- This function calculates the matrices T0 or T1     -%
% %- according to the index (see Vidakovic, 2002)       -%
% %------------------------------------------------------%
% function TT = Tmat(filter, N, index)
%     ind = 1:(N-1);
%     filter2 = zeros(3*N - 5);
%     
%     if index == 0
%         filter2((N-1):(2*N-2)) = sqrt(2)*filter;
%         ft = @(ii, jj) 2*ii - jj - 1;
%     else
%         filter2((N-2):(2*N-3)) = sqrt(2)*filter;
%         ft = @(ii, jj) 2*ii - jj;
%     end
%     
%     ij = bsxfun(ft, ind', ind);
%     ij = ij + 1 - ij(1, end);
%     
%     TT = filter2(ij);
% 
% end
% %---------------------- B. Vidakovic, 2010 --------------------
