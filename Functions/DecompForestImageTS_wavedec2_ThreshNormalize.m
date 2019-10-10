function [ xdec_H, xdec_V, xdec_D, L1 ] = DecompForestImageTS_wavedec2_ThreshNormalize( idx, idy, img_mean, t, J, wname, npts, j1 )
% Description: Here we take logs of the t-th image and of the mean image,
% then we perform a smoothing stage and after that we normalize the t-th
% image. Finally, we take a 2D-DWT, estimate the density function of the
% coefficient's distribution for each detail subband separetely, and store
% the wavelet coefficients of these densities on matrices that are returned
% 
% Input
% idx and idy: indices of the x and y axis that are used in the Forest
% images
% img_mean: mean image of the the Forest images
% t: time of the image to be analyzed
% J: level used in the 2D wavelet decomposition with wavedec2
% wname: name of the wavelet used
% npts: number of points used with the binning estimate of the coefficients densities
% j1: level of the 1D wavelet decomposition of the wavelet densities
% Output
% xdec_H: matrix with wavelet coefficiets of the decomposition of the
% estimated density of horizontal details, for each of the J detail
% subbands
% xdec_V and xdec_D: similar to xdec_H, but for vertical and diagonal details
% L1: number of coefficients on each scale of the wavelet decomposition of
% the coefficient's densities

SatImages = matfile('/media/rodney/Arquivos/Datasets/Abdou_87_Sattelite_Images/xSpatial_ySpatial_VH_VV_Time.mat');

% changing the current boundary to periodic.
dwtmode('per','nodisp');


% Performing the thresholding before normalizing the log-images
mS = log(SatImages.data(idx,idy,1,t) + 1);
[thr,~,~] = ddencmp('cmp','wp',mS);
mG = wdencmp('gbl',mS,'db6',6,thr*2,'s',1);

mS = log(img_mean + 1);
[thr,~,~] = ddencmp('cmp','wp',mS);
mG_mean = wdencmp('gbl',mS,'db6',6,thr*2,'s',1);

% taking the 2D-DWT of the normalization of the thresholded image
[vC,mL] = wavedec2(mG - mG_mean,J,wname);

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

% using the same points of binning for all density estimates
vpts = linspace(-0.5,0.5,npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);

[~,C1,~] = wden(ones(npts,1),'sqtwolog','s','sln',j1,wname);
n_dec = length(C1);

xdec_H = zeros(n_dec,J);
xdec_V = zeros(n_dec,J);
xdec_D = zeros(n_dec,J);

for jj=1:J
    % vectorizing the horizontal detail coefficients stored in cellH
    % binning, which gives a pre-estimate of the square-root of the densities
    bin1 = sqrt(binning(cellH{jj},vpts,wind));
    % wavelet decomposition and denoising of the binned data
    [~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    xdec_H(:,jj) = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

    % performing the same for vertical detail coefficients
    bin1 = sqrt(binning(cellV{jj},vpts,wind));
    [~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    xdec_V(:,jj) = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

    % performing the same for diagonal detail coefficients
    bin1 = sqrt(binning(cellD{jj},vpts,wind));
    [~,C1,L1] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    xdec_D(:,jj) = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));
end


end

