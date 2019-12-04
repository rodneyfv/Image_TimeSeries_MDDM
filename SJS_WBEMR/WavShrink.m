function [yrec, Jest, riskmin] = WavShrink(y, filter, usemad, shrink, wavedomain, signoise)
% WavShrink
% Shrinkage/threshold methods for wavelet coefficients.
%
% Syntax
%   [yrec, Jest, riskmin] = WavShrink(y, filter, usemad, shrink, wavedomain, signoise) 
%   yrec = WavShrink(y, filter, usemad, shrink, wavedomain, signoise) 
%
% Description
% WavShrink applies Efromovich's shrinkage method or Hard/Soft thresholding
%   methods to wavelet coefficient estimates.
%
% Inputs
% y             - Wavelet coefficient estimates. The vector can be the
%                   scaling coefficient estimates or the decomposed wavelet
%                   coefficient estimates.
% filter        - Scaling wavelet filter.
% usemad        - Logical argument. If it is true, the function will
%                   calculate the median absolute deviation (MAD) of the
%                   detail coefficients from the finest level. Otherwise,
%                   the function will calculate the sample standard
%                   deviation of the detail coefficients from the finest
%                   level. Only useful if y corresponds to the scaling
%                   coefficient estimates.
% shrink        - Indicates what method should be used. The possible
%                   methods are 'Hard' for the hard threshold, 'Soft' for
%                   the soft threshold or 'Efromovich' for the Efromovich's
%                   shrinkage method.
% wavedomain    - Logical argument. When it is true, it indicates that the
%                   argument y corresponds to a vector of detail
%                   coefficient estimates (the vector is already in wavelet
%                   domain). When it is false, the method will decompose y
%                   into the wavelet domain and estimate the noise standard
%                   deviation (see usemad argument).
% signoise      - The estimate of the noise standard deviation must be
%                   provided, if the vector y is in wavelet domain.
%
% Outputs
% yrec          - Scaling coefficient estimates reconstructed after
%                   regularization.
% Jest          - Resolution level to be used. This output is useful for
%                   the Efromovich's shrinkage method, and it will
%                   correspond to the resolution level which provides the
%                   minimum empirical risk in our problem.
% riskmin       - The value of the minimum empirical risk.
%
% Details
% The first SINTAX is only useful for the Efromovich's shrinkage method,
%   where it will provide the reconstruced scaling coefficient estimates
%   (yrec), the selected resolution level that minimizes the empirical risk
%   (Jest) and the minimum empirical risk (riskmin). The second sintax can
%   be used for all the three methodologies of regularization, and it
%   returns only the reconstruced coefficients.
%
% Example
% % See its use in the Application.m file

% By Michel H. Montoril at Georgia Institute of Technology, Sep 2014.

    n = length(y);
    
    if wavedomain
        wc = y;
        sig = signoise;
    else
        [wc, sig] = wd(y, 0, filter, usemad);
    end
    
    if strcmp(shrink, 'Hard')
        Jest = log2(n);
        riskmin = NaN;
        
        lambda = sig*sqrt(2*log(n));
        wct = wthresh(wc,'h',lambda);
        yrec = wr(wct, 0, filter);
        
    elseif strcmp(shrink, 'Soft')
        Jest = log2(n);
        riskmin = NaN;
        
        lambda = sig*sqrt(2*log(n));
        wct = wthresh(wc,'s',lambda);
        yrec = wr(wct, 0, filter);
        
    elseif strcmp(shrink, 'Efromovich')
        J = log2(n)-log2(2-wavedomain);
        THETA = max(wc(1:2^J).^2-sig^2, 0);
        
        s1 = (1:n/(2-wavedomain))*sig^2;
        s2 = sum(THETA) - cumsum(THETA);
        emprisk = s1 + s2;
        wemprisk = NaN(1, 2^J);
        wemprisk(1:2) = emprisk(1:2);
        
        if J > 1
            for kk= 2:J
                wemprisk(kk+1) = mean(emprisk((2^(kk-1)+1):(2^kk)));
            end
        end
        
        [riskmin, Jest] = min(wemprisk);
    
        Jest = Jest - 1;
        ind = 1:2^Jest;
        wc = wc(ind);
        THETA = THETA(ind);
        cTHETA = THETA ./ (THETA + sig^2);
    
        shrinkwc = cTHETA .* wc;
        yrec = wr(shrinkwc, 0, filter);
        
    end
    
end