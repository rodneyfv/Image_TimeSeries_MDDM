function [dist,distE]=Hellinger_WT_interp(im1,im2,wname,J)
% This function computes the total Hellinger distance of two images im1 and
% im2. The density function of the wavelet coefficents in each sub-matrix
% of the decomposition of these images, is estimated non-parametricaly
% using the wavelet density estimator proposed by Pinheiro and Vidakovic
% (1997). Here, the wavelet functions employed are computed based on pre-
% computed values of the scaling and wavelet functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decomposition de l'approximation swpt par une swt
im1 = im1 - mean2(im1);
im2 = im2 - mean2(im2);
[SCA1,SCH1,SCV1,SCD1]=swt2(im1,J,wname);
[SCA2,SCH2,SCV2,SCD2]=swt2(im2,J,wname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = load('sym8_phi_values.mat');
mValues_phi = tmp.mValues;
tmp = load('sym8_psi_values.mat');
mValues_psi = tmp.mValues;

N = length(wfilters(wname))/2;
j0 = ceil(log2(2*N-1));
j1 = 5;

% vector of points equally spaced and in the same interval that
% contains both values of img1 and img2
vpts = linspace(0,1,500);

dist=0;
% Traitement des details swt
for j=1:J
    img1=SCV1(:,:,j);
    img2=SCV2(:,:,j);
    [v_alpha, v_beta] = dwt_01_interp( img1(:), 2, wname, j0, j1, mValues_phi, mValues_psi);
    sqrtdens1 = sqrtdens_dwt_01_interp( vpts, v_alpha, v_beta, wname, j0, j1, mValues_phi, mValues_psi);
    [v_alpha, v_beta] = dwt_01_interp( img2(:), 2, wname, j0, j1, mValues_phi, mValues_psi);
    sqrtdens2 = sqrtdens_dwt_01_interp( vpts, v_alpha, v_beta, wname, j0, j1, mValues_phi, mValues_psi);
    % computing Hellinger distance
    dist = dist + sqrt(0.5*sum((sqrtdens1 - sqrtdens2).^2)*(vpts(2)-vpts(1)));
    % the same process is repeated for the other submatrices of the SWT
    % decomposition
    img1=SCH1(:,:,j);
    img2=SCH2(:,:,j);
    [v_alpha, v_beta] = dwt_01_interp( img1(:), 2, wname, j0, j1, mValues_phi, mValues_psi);
    sqrtdens1 = sqrtdens_dwt_01_interp( vpts, v_alpha, v_beta, wname, j0, j1, mValues_phi, mValues_psi);
    [v_alpha, v_beta] = dwt_01_interp( img2(:), 2, wname, j0, j1, mValues_phi, mValues_psi);
    sqrtdens2 = sqrtdens_dwt_01_interp( vpts, v_alpha, v_beta, wname, j0, j1, mValues_phi, mValues_psi);
    dist = dist + sqrt(0.5*sum((sqrtdens1 - sqrtdens2).^2)*(vpts(2)-vpts(1)));
    img1=SCD1(:,:,j);
    img2=SCD2(:,:,j);
    [v_alpha, v_beta] = dwt_01_interp( img1(:), 2, wname, j0, j1, mValues_phi, mValues_psi);
    sqrtdens1 = sqrtdens_dwt_01_interp( vpts, v_alpha, v_beta, wname, j0, j1, mValues_phi, mValues_psi);
    [v_alpha, v_beta] = dwt_01_interp( img2(:), 2, wname, j0, j1, mValues_phi, mValues_psi);
    sqrtdens2 = sqrtdens_dwt_01_interp( vpts, v_alpha, v_beta, wname, j0, j1, mValues_phi, mValues_psi);
    dist = dist + sqrt(0.5*sum((sqrtdens1 - sqrtdens2).^2)*(vpts(2)-vpts(1)));
end
% performing the same process as before for the approximation coefficients
% of the SWT decomposition
img1=SCA1(:,:,j);
img2=SCA2(:,:,j);
[v_alpha, v_beta] = dwt_01_interp( img1(:), 2, wname, j0, j1, mValues_phi, mValues_psi);
sqrtdens1 = sqrtdens_dwt_01_interp( vpts, v_alpha, v_beta, wname, j0, j1, mValues_phi, mValues_psi);
[v_alpha, v_beta] = dwt_01_interp( img2(:), 2, wname, j0, j1, mValues_phi, mValues_psi);
sqrtdens2 = sqrtdens_dwt_01_interp( vpts, v_alpha, v_beta, wname, j0, j1, mValues_phi, mValues_psi);
distE = sqrt(0.5*sum((sqrtdens1 - sqrtdens2).^2)*(vpts(2)-vpts(1)));
        
end
