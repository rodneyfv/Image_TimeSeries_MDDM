function [ xdec_A, xdec_H, xdec_V, xdec_D, L1 ] = DecompForestImageTS_swt2( idx, idy, t, J, wname, npts, j1 )
% With this function we access the Forest data, get the image corresponding
% to time t and at the pixels idx and idy, then perform the MDDM process
% using the swt2 function for the 2D wavelet decomposition of the image. J
% is the level of swt2 transform, wname is the wavelet type used, npts is
% the number of points used in the estimated density with the binning and
% j1 is the level used in the thresholding of the estimated density.

SatImages = matfile('/media/rodney/Arquivos/Datasets/Abdou_87_Sattelite_Images/xSpatial_ySpatial_VH_VV_Time.mat');

img = SatImages.data(idx,idy,1,t);

mS = log(img + 1);
% smoothing the image using 2D wavelet denoising
[thr,sorh,keepapp] = ddencmp('cmp','wp',mS);
mG = wdencmp('gbl',mS,'db6',6,thr*2,'s',1);
[SCA1,SCH1,SCV1,SCD1]=swt2(mG - mean2(mG),J,wname);

% using the same points of binning for all density estimates
vpts = linspace(-0.5,0.5,npts);

img = SCA1;
% points where the binning will be performed
%vpts = linspace(min(min(img)),max(max(img)),npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);
% binning, which gives a pre-estimate of the square-root of the densities
bin1 = sqrt(binning(img(:),vpts,wind));
% wavelet decomposition and denoising of the binned data
[~,C1,L1] = wden(bin1,'sqtwolog','s','sln',j1,wname);
xdec_A = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

img = SCH1;
% points where the binning will be performed
%vpts = linspace(min(min(img)),max(max(img)),npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);
% binning, which gives a pre-estimate of the square-root of the densities
bin1 = sqrt(binning(img(:),vpts,wind));
% wavelet decomposition and denoising of the binned data
[~,C1,L1] = wden(bin1,'sqtwolog','s','sln',j1,wname);
xdec_H = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

img = SCV1;
% points where the binning will be performed
%vpts = linspace(min(min(img)),max(max(img)),npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);
% binning, which gives a pre-estimate of the square-root of the densities
bin1 = sqrt(binning(img(:),vpts,wind));
% wavelet decomposition and denoising of the binned data
[~,C1,L1] = wden(bin1,'sqtwolog','s','sln',j1,wname);
xdec_V = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

img = SCD1;
% points where the binning will be performed
%vpts = linspace(min(min(img)),max(max(img)),npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);
% binning, which gives a pre-estimate of the square-root of the densities
bin1 = sqrt(binning(img(:),vpts,wind));
% wavelet decomposition and denoising of the binned data
[~,C1,L1] = wden(bin1,'sqtwolog','s','sln',j1,wname);
xdec_D = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

end

