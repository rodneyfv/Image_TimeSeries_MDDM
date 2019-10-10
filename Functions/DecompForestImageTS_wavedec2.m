function [ xdec_H, xdec_V, xdec_D, L1 ] = DecompForestImageTS_wavedec2( idx, idy, img_mean, t, J, J_MDDM, wname, npts, j1 )
% With this function we access the Forest data, get the image corresponding
% to time t and at the pixels idx and idy, then perform the MDDM process
% using the wavedec2 function for the 2D wavelet decomposition of the
% image minus the mean image, which is stored in img_mean. J is the level 
% of wavedec2 transform, J_MDDM is the maximum level
% of detail coefficients we analyze when computing the density function of 
% the wavelet coefficients, wname is the wavelet type used, npts is
% the number of points used in the estimated density with the binning and
% j1 is the level used in the thresholding of the estimated density.

SatImages = matfile('/media/rodney/Arquivos/Datasets/Abdou_87_Sattelite_Images/xSpatial_ySpatial_VH_VV_Time.mat');

img = SatImages.data(idx,idy,1,t) - img_mean;

mS = log(img + 1);
% smoothing the image using 2D wavelet denoising
[thr,~,~] = ddencmp('cmp','wp',mS);
mG = wdencmp('gbl',mS,'db6',6,thr*2,'s',1);
[vC,mL] = wavedec2(mG - mean2(mG),J,wname);

% getting the matrices with detail coefficients up to level J_MDDM and
% storing them as a vector in a cell
cellH = cell(J_MDDM);
cellV = cell(J_MDDM);
cellD = cell(J_MDDM);
for jj=1:J_MDDM
    [tmpH,tmpV,tmpD] = detcoef2('all',vC,mL,jj);
    cellH{jj} = tmpH(:);
    cellV{jj} = tmpV(:);
    cellD{jj} = tmpD(:);
end

% using the same points of binning for all density estimates
vpts = linspace(-0.5,0.5,npts);

% vectorizing the horizontal detail coefficients stored in cellH
img = cell2mat(cellH);
% points where the binning will be performed
%vpts = linspace(min(img),max(img),npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);
% binning, which gives a pre-estimate of the square-root of the densities
bin1 = sqrt(binning(img,vpts,wind));
% wavelet decomposition and denoising of the binned data
[~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
xdec_H = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

% performing the same for vertical detail coefficients
img = cell2mat(cellV);
%vpts = linspace(min(img),max(img),npts);
wind = (vpts(npts) - vpts(1))/sqrt(npts);
bin1 = sqrt(binning(img,vpts,wind));
[~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
xdec_V = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

% performing the same for diagonal detail coefficients
img = cell2mat(cellD);
%vpts = linspace(min(img),max(img),npts);
wind = (vpts(npts) - vpts(1))/sqrt(npts);
bin1 = sqrt(binning(img,vpts,wind));
[~,C1,L1] = wden(bin1,'sqtwolog','s','sln',j1,wname);
xdec_D = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

end

