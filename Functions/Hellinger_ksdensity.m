function [dist,distE]=Hellinger_ksdensity(im1,im2,wname,J)
% This function computes the total Hellinger distance of two images im1 and
% im2. The density function of the wavelet coefficents in each sub-matrix
% of the decomposition of these images, is estimated non-parametricaly
% using a kernel estimator.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decomposition de l'approximation swpt par une swt
im1 = im1 - mean2(im1);
im2 = im2 - mean2(im2);
[SCA1,SCH1,SCV1,SCD1]=swt2(im1,J,wname);
[SCA2,SCH2,SCV2,SCD2]=swt2(im2,J,wname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist=0;
% Traitement des details swt
for j=1:J
    % vertical details
    img1=SCV1(:,:,j);
    img2=SCV2(:,:,j);
    % vector of points equally spaced and in the same interval that
    % contains both values of img1 and img2
    vpts = linspace(min(min(img1(:)),min(img2(:))),max(max(img1(:)),max(img2(:))),1000);
    % density estimate of the values on both vectorized image matrices
    densimg1 = ksdensity(img1(:),vpts);
    densimg2 = ksdensity(img2(:),vpts);
    % computing Hellinger distance
    dist = dist + 1 - sum(sqrt(densimg1.*densimg2))*(vpts(2)-vpts(1));
    % the same process is repeated for the other submatrices of the SWT
    % decomposition
    % horizontal details
    img1=SCH1(:,:,j);
    img2=SCH2(:,:,j);
    vpts = linspace(min(min(img1(:)),min(img2(:))),max(max(img1(:)),max(img2(:))),1000);
    densimg1 = ksdensity(img1(:),vpts);
    densimg2 = ksdensity(img2(:),vpts);
    dist = dist + 1 - sum(sqrt(densimg1.*densimg2))*(vpts(2)-vpts(1));
    % diagonal details
    img1=SCD1(:,:,j);
    img2=SCD2(:,:,j);
    vpts = linspace(min(min(img1(:)),min(img2(:))),max(max(img1(:)),max(img2(:))),1000);
    densimg1 = ksdensity(img1(:),vpts);
    densimg2 = ksdensity(img2(:),vpts);
    dist = dist + 1 - sum(sqrt(densimg1.*densimg2))*(vpts(2)-vpts(1));
end
% performing the same process as before for the approximation coefficients
% of the SWT decomposition
img1=SCA1(:,:,J);
img2=SCA2(:,:,J);
vpts = linspace(min(min(img1(:)),min(img2(:))),max(max(img1(:)),max(img2(:))),1000);
densimg1 = ksdensity(img1(:),vpts);
densimg2 = ksdensity(img2(:),vpts);
distE = 1 - sum(sqrt(densimg1.*densimg2))*(vpts(2)-vpts(1));
        
end