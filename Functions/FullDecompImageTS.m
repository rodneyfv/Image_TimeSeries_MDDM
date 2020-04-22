function [ x_A, x_H, x_V, x_D ] = FullDecompImageTS( img_obs, J, wname, npts, j1 )
%
% Description
% 
% We apply the following steps, in this order: a wavelet denoising of the
% image img_obs if performed; a 2D-DWT decomposition is performed with J
% levels; for each kind of detail (approximation, horizontal, vertical,
% diagonal) coefficient, a density function of the coefficients is
% estimated through binning; a 1D-DWT wavelet denoising of the values
% of this density is performed and returned in the function output.
% Pay attention that values where binned is performed is different for
% approximation coefficients than the other detail coefficients.
% 
% Inputs
% img_obs   - image to be analyzed (logs are not taken inside this function)
% J         - number of decomposition levels taken in the 2D-DWT of the image
% wname     - wavelet used
% npts      - number of points used in the binning when estimating the density
%           function of the 2D-DWT coefficients
% j1        - number of decomposition levels taken in the 1D-DWT of the binning
%           estimate of the density
% 
% Outputs
% x_A    - denoised binning density estimate for the coefficients of the 
%       approximation level of the image 2D-DWT decomposition
% x_H    - ... of the horizontal details of the image 2D-DWT decomposition
% x_V    - ... of the vertical details of the image 2D-DWT decomposition
% x_D    - ... of the diagonal details of the image 2D-DWT decomposition

% changing the current boundary to periodic.
dwtmode('per','nodisp');

% smoothing the image using 2D wavelet denoising
[thr,~,~] = ddencmp('cmp','wp',img_obs);
mG = wdencmp('gbl',img_obs,'db6',6,thr*2,'s',1);

% taking the 2D-DWT of the thresholded normalized image
[vC,mL] = wavedec2(mG,J,wname);


% getting the matrices with detail coefficients up to level J_MDDM and
% storing them as a vector in a cell
cellH = cell(J);
cellV = cell(J);
cellD = cell(J);
for jj=1:J
    [tmpH,tmpV,tmpD] = detcoef2('all',vC,mL,jj);
    cellH{jj} = tmpH(:);
    cellV{jj} = tmpV(:);
    cellD{jj} = tmpD(:);
end
% vector with approximation coefficients
cellA = vC(1:prod(mL(1,:)));

% using the same points of binning for all density estimates
vpts = linspace(-0.5,0.5,npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);

x_H = zeros(npts,J);
x_V = zeros(npts,J);
x_D = zeros(npts,J);

for jj=1:J
    % vectorizing the horizontal detail coefficients stored in cellH
    % binning, which gives a pre-estimate of the square-root of the densities
    bin1 = sqrt(binning(cellH{jj},vpts,wind));
    % wavelet decomposition and denoising of the binned data
    [x1,~,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    x_H(:,jj) = x1/sqrt(sum(x1.^2)*(vpts(2)-vpts(1)));

    % performing the same for vertical detail coefficients
    bin1 = sqrt(binning(cellV{jj},vpts,wind));
    [x1,~,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    x_V(:,jj) = x1/sqrt(sum(x1.^2)*(vpts(2)-vpts(1)));

    % performing the same for diagonal detail coefficients
    bin1 = sqrt(binning(cellD{jj},vpts,wind));
    [x1,~,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    x_D(:,jj) = x1/sqrt(sum(x1.^2)*(vpts(2)-vpts(1)));    
end
% using different points of binning for approx. coeff. density
vpts = linspace(-5,40,npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);
% performing the same for approximation coefficients
bin1 = sqrt(binning(cellA,vpts,wind));
[x1,~,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
x_A = x1/sqrt(sum(x1.^2)*(vpts(2)-vpts(1)));


end

