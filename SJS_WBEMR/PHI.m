function [PhiMat, rangek] = PHI(x, J, filter, prec, periodized)
% PHI
% Daubechies-Lagarias approximation of scaling functions for a data set.
%
% Syntax
%   [PhiMat, rangek] = PHI(x, J, filter, prec, periodized)
%   PhiMat = PHI(x, J, filter, prec, periodized)
%
% Description
% PHI calculates the scaling functions corresponding to an orthogonal MRA
%   by Daubechies-Lagarias algorithm. Basically, the function calculates $
%   phi_{Jk}(x_i) $, where $ x_i $ corresponds to the i-th datum of a data
%   set x.
%
% Inputs
% x           - Data set to calculate the scaling function.
% J           - Resolution level.
% filter      - Scaling wavelet filter.
% prec        - Precision approximation measured by the number of
%                 Daubechies-Lagarias steps.
% preriodized - Logial argument. Should the scaling functions be
%                 periodized? If it is true, the scaling functions will be
%                 periodized. We recommend only to use true.
%
% Outputs
% PhiMat      - Matrix of scaling functions. Each line corresponds
%                 to an element of x, and each column corresponds to an
%                 element $ k $ of the basis $ \{ \phi_{Jk} \}_k $. The
%                 number of columns will depend on the wavelet filter to be
%                 used and if the periodized argument is true of false.
% rangek      - Assuming that the support of the scaling function $ \phi $
%                 is in $ [0, N-1] $, where $ N $ is the
%                 length of the scaling wavelet filter, this output
%                 gives the range of variation of $ k $.
%
% Example
% % In this example we consider:
% % A random data set in the interval (0, 1) called data
% data = rand(5, 1);
% % A Daubechies Least Asymmetric 10-tap wavelet filter called wfilter
% wfilter = [0.027333068345164, 0.029519490926073,-0.039134249302583,...
%            0.199397533976996, 0.723407690403808, 0.633978963456949,...
%            0.016602105764424,-0.175328089908107,-0.021101834024930,...
%            0.019538882735387];
% % A resolution level called J
% J = 3;
% % A number of Daubechies-Lagarias steps prec
% prec = 30;
% % The resulting matrix of periodized scaling functions (matPHI) and
% % the range of variarion of k (krange) before periodization are
% % calculated by
% [matPHI, krange] = PHI(data, J, wfilter, 30, true);
%

% By Michel H. Montoril at Georgia Institute of Technology, Sep 2014.

   %#codegen
   coder.inline('never')

   n = length(x);
   N = length(filter);
   p = 2^J;
   y = p*x;
   kmin = ceil(0*min(y) - N + 1);
   kmax = floor(max(y) - 1e-9);
   rangek = [kmin kmax];
   PhiMat = zeros(n, kmax - kmin + 1);
   for i = 1:n
       k1 = ceil(y(i) - N + 1);
       PhiMat(i, (k1 - kmin + N - 1):-1:(k1 - kmin + 1)) = sqrt(p)*Phi(y(i), filter, N, prec);
   end
   
   if periodized == true
       PhiMat = sum(reshape([PhiMat zeros(n, p - mod(size(PhiMat, 2), p))], n, p, []), 3);
       PhiMat = PhiMat(:, mod(1:p, p) + 1);
       %PhiMat = PhiMat(:, mod(kmin+N:(kmin+N + p - 1), p) + 1);
   end
end

function phi = Phi(x, filter, N, prec)
    int = floor(x);
    dec = x - int;
    dy = dec2bin(dec, prec);
    T0 = Tmat(filter, N, 0);
    T1 = Tmat(filter, N, 1);
    prod = eye(N - 1);
    for i = 1:prec
        if dy(i)==1
            prod = prod*T1;
        else
            prod = prod*T0;
        end
    end
    phi = mean(prod, 2);
end

function a = dec2bin(x, n)
    a = zeros(1, n);
    y = x;
    for i = 1:n
         y = 2*y;
        a(i) = floor(y);
        y = y - floor(y);
    end
end

function TT = Tmat(filter, N, index)
    ind = 1:(N-1);
    filter2 = zeros(3*N - 5);
    
    if index == 0
        filter2((N-1):(2*N-2)) = sqrt(2)*filter;
        ij = bsxfun(@ft0, ind', ind);
    else
        filter2((N-2):(2*N-3)) = sqrt(2)*filter;
        ij = bsxfun(@ft1, ind', ind);
    end
    
    ij = ij + 1 - ij(1, end);
    TT = filter2(ij);

end

function val0 = ft0(ii, jj)
    val0 = 2*ii - jj - 1;
end

function val1 = ft1(ii, jj)
    val1 = 2*ii - jj;
end