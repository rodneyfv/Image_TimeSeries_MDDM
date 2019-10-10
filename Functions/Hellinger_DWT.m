function [dist,distE] = Hellinger_DWT(im1,im2,wname,J)
% This function computes the total Hellinger distance of two images im1 and
% im2. A binning method is applied to the data, and an approximation of the
% density function is computed at equispaced points. The values obtained
% are decomposed with the DWT, whose wavelet coefficients generated are
% used to compute the Hellinger distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decomposition de l'approximation swpt par une swt
im1 = im1 - mean2(im1);
im2 = im2 - mean2(im2);
[SCA1,SCH1,SCV1,SCD1]=swt2(im1,J,wname);
[SCA2,SCH2,SCV2,SCD2]=swt2(im2,J,wname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of equispaced points to evaluate the binned data
npts = min(512,2^floor(log2(numel(im1))-4));
% level of the 1D-DWT used for the equispaced values
j1 = floor(log2(npts) - log2(log(npts))); % See Eq. (18) of PV1997

dist=0;
% Traitement des details swt
for jj=1:J
    img1=SCV1(:,:,jj);
    img2=SCV2(:,:,jj);
    % points where the binning will be performed
    vpts = linspace(min([min(img1),min(img2)]),max([max(img1),max(img2)]),npts);
    % window length used in the binning
    wind = (vpts(npts) - vpts(1))/sqrt(npts);
    % binning, which gives a pre-estimate of the square-root of the densities
    bin1 = sqrt(binning(img1(:),vpts,wind));
    bin2 = sqrt(binning(img2(:),vpts,wind));
    % wavelet decomposition and denoising of the binned data
    [~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    [~,C2,~] = wden(bin2,'sqtwolog','s','sln',j1,wname);
    % updating with the Hellinger divergence
    dist = dist + norm(C1/norm(C1) - C2/norm(C2))/sqrt(2);    
    
    % the same process is repeated for the other submatrices of the 2D-DWT
    % decomposition
    img1=SCH1(:,:,jj);
    img2=SCH2(:,:,jj);
    vpts = linspace(min([min(img1),min(img2)]),max([max(img1),max(img2)]),npts);
    wind = (vpts(npts) - vpts(1))/sqrt(npts);
    bin1 = sqrt(binning(img1(:),vpts,wind));
    bin2 = sqrt(binning(img2(:),vpts,wind));
    [~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    [~,C2,~] = wden(bin2,'sqtwolog','s','sln',j1,wname);
    dist = dist + norm(C1/norm(C1) - C2/norm(C2))/sqrt(2);
    
    img1=SCD1(:,:,jj);
    img2=SCD2(:,:,jj);
    vpts = linspace(min([min(img1),min(img2)]),max([max(img1),max(img2)]),npts);
    wind = (vpts(npts) - vpts(1))/sqrt(npts);
    bin1 = sqrt(binning(img1(:),vpts,wind));
    bin2 = sqrt(binning(img2(:),vpts,wind));
    [~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    [~,C2,~] = wden(bin2,'sqtwolog','s','sln',j1,wname);
    dist = dist + norm(C1/norm(C1) - C2/norm(C2))/sqrt(2);
end
% performing the same process as before for the approximation coefficients
% of the 2D-DWT decomposition
img1=SCA1(:,:,jj);
img2=SCA2(:,:,jj);
vpts = linspace(min([min(img1),min(img2)]),max([max(img1),max(img2)]),npts);
wind = (vpts(npts) - vpts(1))/sqrt(npts);
bin1 = sqrt(binning(img1(:),vpts,wind));
bin2 = sqrt(binning(img2(:),vpts,wind));
[~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
[~,C2,~] = wden(bin2,'sqtwolog','s','sln',j1,wname);
distE = norm(C1/norm(C1) - C2/norm(C2))/sqrt(2);

        
end
