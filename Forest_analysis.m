% 1 - Here we test the idea of estimating the dimension of the subspace
% generating the density time series of wavelet coefficients, then using
% the densities estimated with the eigenfunctions and the dimension chosen
% to obtain residuals, which are employed to obtain a bootstrap
% distribution of the Hellinger distance under H0.
% 2 - Here we make the same analysis as before, but also using the
% approximation coefficients.

close all
clear all

addpath(genpath('Functions/'))

SatImages = matfile('/media/rodney/Arquivos/Datasets/Abdou_87_Sattelite_Images/xSpatial_ySpatial_VH_VV_Time.mat');
% SatImages.data is a matrix of dimension 4000x3000x2x87

%ds = 16;
%idx = 1:ds:2048;
%idy = 501:ds:2548;
idx = 1:1024;
idy = 501:(500 + 1024);

%img = SatImages.data(idx,idy,1,1);

wname='sym8';
%wname='db4';

% saving the current boundary handling and then changing it to periodic.
% The default mode was giving different norms of vC for different
% decomposition levels J
origMode = dwtmode('status','nodisplay');
dwtmode('per');

% number of equispaced points to evaluate the binned data
%npts = 1024;
npts = 512;
% level of the 1D-DWT used for the equispaced values
j1 = floor(log2(npts) - log2(log(npts))); % See Eq. (18) of PV1997
[~,tmp,~] = wden(ones(npts,1),'sqtwolog','s','sln',j1,wname);
n_dec = length(tmp);

% level chosen using the scalogram
J = 6;

% sample size of the time series
n = 87;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we test the idea of estimating the dimension of the subspace
% generating the density time series of wavelet coefficients, then using
% the densities estimated with the eigenfunctions and the dimension chosen
% to obtain residuals, which are employed to obtain a bootstrap
% distribution of the Hellinger distance under H0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computing the mean image
img_mean = zeros(length(idx), length(idy));
tic
for t=1:5:86
    % we compute directly the sum of five images each time
    if t<=85
        img = SatImages.data(idx,idy,1,t:(t+4));
        img_mean = img_mean + sum(img,4)./n;
    else
        img = SatImages.data(idx,idy,1,t:n);
        img_mean = img_mean + sum(img,4)./n;        
    end
end
toc


% arrays with wavelet coefficients of densities for each of the J details
Xdec_H = zeros(n_dec,J,n);
Xdec_V = zeros(n_dec,J,n);
Xdec_D = zeros(n_dec,J,n);
parfor jj=1:2 % initiating the parfor
    tmp = 0;
end
parfor t=1:n
    [Xdec_H(:,:,t), Xdec_V(:,:,t), Xdec_D(:,:,t), ~] = DecompForestImageTS_wavedec2_NormalizeThresh( idx, idy, img_mean, t, J, wname, npts, j1 );
end


% estimating the dimension of the density time series for each sub-region
% of the 2D decomposition of the images
p = 2;  % maximum lag used
Nboot = 301;  % number of bootstrap replications
alpha = .05;  % significance level of the bootstrap test

[dim_H,vPvalues_H,~] = Estim_Dim_Pval( reshape(Xdec_H(:,1,:),[npts n]), p, Nboot, alpha );
dim_H
vPvalues_H
[dim_V,vPvalues_V,~] = Estim_Dim_Pval( reshape(Xdec_V(:,1,:),[npts n]), p, Nboot, alpha );
dim_V
vPvalues_V
[dim_D,vPvalues_D,~] = Estim_Dim_Pval( reshape(Xdec_D(:,1,:),[npts n]), p, Nboot, alpha );
dim_D
vPvalues_D


% estimated dimension
d = 1;

% matrices for the wavelet coefficients of the densities corresponding to the
% estimated dimension
Xdec_hat_H = zeros(n_dec,J,n);
Xdec_hat_V = zeros(n_dec,J,n);
Xdec_hat_D = zeros(n_dec,J,n);
% matrices of residuals obtained from the difference between the observed
% densities and the estimated densities
Eps_hat_H = zeros(J,n);
Eps_hat_V = zeros(J,n);
Eps_hat_D = zeros(J,n);

for jj=1:J
    % estimated density
    Xdec_hat_H(:,jj,:) = estimated_functions_Dim( reshape(Xdec_H(:,jj,:),[npts n]), p, d );
    Xdec_hat_V(:,jj,:) = estimated_functions_Dim( reshape(Xdec_V(:,jj,:),[npts n]), p, d );
    Xdec_hat_D(:,jj,:) = estimated_functions_Dim( reshape(Xdec_D(:,jj,:),[npts n]), p, d );
    % residual obtained
    for t=1:n
        Eps_hat_H(jj,t) = norm(Xdec_H(:,jj,t) - Xdec_hat_H(:,jj,t));
        Eps_hat_V(jj,t) = norm(Xdec_V(:,jj,t) - Xdec_hat_V(:,jj,t));
        Eps_hat_D(jj,t) = norm(Xdec_D(:,jj,t) - Xdec_hat_D(:,jj,t));
    end
end


% computing the MDDMs
MDDM_V = zeros(n,n,J);
MDDM_H = zeros(n,n,J);
MDDM_D = zeros(n,n,J);
for jj=1:J
    for k=1:(n-1)
        for m=(k+1):n
            MDDM_H(k,m,jj) = norm(Xdec_H(:,jj,k) - Xdec_H(:,jj,m))/sqrt(2);
            MDDM_H(m,k,jj) = MDDM_H(k,m,jj);
            MDDM_V(k,m,jj) = norm(Xdec_V(:,jj,k) - Xdec_V(:,jj,m))/sqrt(2);
            MDDM_V(m,k,jj) = MDDM_V(k,m,jj);
            MDDM_D(k,m,jj) = norm(Xdec_D(:,jj,k) - Xdec_D(:,jj,m))/sqrt(2);
            MDDM_D(m,k,jj) = MDDM_D(k,m,jj);
        end
    end    
end

% here we use the residuals obtained to generate bootstrap samples of the
% Hellinger distance under H0
vMDDM_boot_H = randsample(Eps_hat_H(5,:),1000,true);
tmp = MDDM_H(:,:,5)>quantile(vMDDM_boot_H,.999);
imagesc(tmp)

vMDDM_boot_V = randsample(Eps_hat_V(5,:),1000,true);
tmp = MDDM_V(:,:,5)>quantile(vMDDM_boot_V,.999);
imagesc(tmp)

vMDDM_boot_H = randsample(Eps_hat_H(5,:),1000,true);
tmp = MDDM_H(:,:,5)>quantile(vMDDM_boot_H,.999);
imagesc(tmp)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we make the same analysis as before, but also using the approximation
% coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computing the mean image
img_mean = zeros(length(idx), length(idy));
for t=1:5:86
    % we compute directly the sum of five images each time
    if t<=85
        img = SatImages.data(idx,idy,1,t:(t+4));
        img_mean = img_mean + sum(img,4)./n;
    else
        img = SatImages.data(idx,idy,1,t:n);
        img_mean = img_mean + sum(img,4)./n;        
    end
end


% arrays with wavelet coefficients of densities for each of the J details
Xdec_A = zeros(n_dec,n);
Xdec_H = zeros(n_dec,J,n);
Xdec_V = zeros(n_dec,J,n);
Xdec_D = zeros(n_dec,J,n);
parfor jj=1:2 % initiating the parfor
    tmp = 0;
end
parfor t=1:n
    %SatImages = matfile('/media/rodney/Arquivos/Datasets/Abdou_87_Sattelite_Images/xSpatial_ySpatial_VH_VV_Time.mat');
    % changing the current boundary to periodic.
    dwtmode('per','nodisp');
    % Normalizing the log-images before performing the thresholding
    mS = log(SatImages.data(idx,idy,1,t) + 1) - log(img_mean + 1);
    [Xdec_A(:,t), Xdec_H(:,:,t), Xdec_V(:,:,t), Xdec_D(:,:,t), ~] = FullDecompImageTS_wavedec2( mS, J, wname, npts, j1 );
end


% computing the MDDMs
MDDM_A = zeros(n,n);
MDDM_V = zeros(n,n,J);
MDDM_H = zeros(n,n,J);
MDDM_D = zeros(n,n,J);
for k=1:(n-1)
    for m=(k+1):n
        MDDM_A(k,m) = norm(Xdec_A(:,k) - Xdec_A(:,m))/sqrt(2);
        MDDM_A(m,k) = MDDM_A(k,m);
        for jj=1:J
            MDDM_H(k,m,jj) = norm(Xdec_H(:,jj,k) - Xdec_H(:,jj,m))/sqrt(2);
            MDDM_H(m,k,jj) = MDDM_H(k,m,jj);
            MDDM_V(k,m,jj) = norm(Xdec_V(:,jj,k) - Xdec_V(:,jj,m))/sqrt(2);
            MDDM_V(m,k,jj) = MDDM_V(k,m,jj);
            MDDM_D(k,m,jj) = norm(Xdec_D(:,jj,k) - Xdec_D(:,jj,m))/sqrt(2);
            MDDM_D(m,k,jj) = MDDM_D(k,m,jj);
        end
    end    
end

imagesc(MDDM_A)
imagesc(sum(MDDM_H + MDDM_V + MDDM_D,3))


% estimating the dimension of the density time series for each sub-region
% of the 2D decomposition of the images
p = 2;  % maximum lag used
Nboot = 200;  % number of bootstrap replications
alpha = .05;  % significance level of the bootstrap test

dim_A = 0;
vdim_H = ones(J,1);
vdim_V = ones(J,1);
vdim_D = ones(J,1);
[dim_A,~,~] = Estim_Dim_Pval( reshape(Xdec_A,[npts n]), p, Nboot, alpha );
parfor jj=1:J
    [vdim_H(jj),~,~] = Estim_Dim_Pval( reshape(Xdec_H(:,jj,:),[npts n]), p, Nboot, alpha );
    [vdim_V(jj),~,~] = Estim_Dim_Pval( reshape(Xdec_V(:,jj,:),[npts n]), p, Nboot, alpha );
    [vdim_D(jj),~,~] = Estim_Dim_Pval( reshape(Xdec_D(:,jj,:),[npts n]), p, Nboot, alpha );
end

% matrices for the wavelet coefficients of the densities corresponding to the
% estimated dimension
Xdec_hat_H = zeros(n_dec,J,n);
Xdec_hat_V = zeros(n_dec,J,n);
Xdec_hat_D = zeros(n_dec,J,n);
% matrices of estimated Hellinger distances under H0, obtained from the 
% difference between the observed densities and the estimated densities
Hel_hat_A = zeros(1,n);
Hel_hat_H = zeros(J,n);
Hel_hat_V = zeros(J,n);
Hel_hat_D = zeros(J,n);

for jj=1:J
    % estimated density
    Xdec_hat_H(:,jj,:) = estimated_functions_Dim( reshape(Xdec_H(:,jj,:),[npts n]), p, vdim_H(jj) );
    Xdec_hat_V(:,jj,:) = estimated_functions_Dim( reshape(Xdec_V(:,jj,:),[npts n]), p, vdim_V(jj) );
    Xdec_hat_D(:,jj,:) = estimated_functions_Dim( reshape(Xdec_D(:,jj,:),[npts n]), p, vdim_D(jj) );
    % residual obtained
    for t=1:n
        Hel_hat_H(jj,t) = norm(Xdec_H(:,jj,t) - Xdec_hat_H(:,jj,t));
        Hel_hat_V(jj,t) = norm(Xdec_V(:,jj,t) - Xdec_hat_V(:,jj,t));
        Hel_hat_D(jj,t) = norm(Xdec_D(:,jj,t) - Xdec_hat_D(:,jj,t));
    end
end
Xdec_hat_A = estimated_functions_Dim( reshape(Xdec_A,[npts n]), p, dim_A );
for t=1:n; Hel_hat_A(t) = norm(Xdec_A(:,t) - Xdec_hat_A(:,t)); end


% here we use the residuals obtained to generate bootstrap samples of the
% Hellinger distance under H0
vMDDM_boot = randsample(sum(Hel_hat_H,1) + sum(Hel_hat_V,1) + sum(Hel_hat_D,1) + Hel_hat_A,1000,true);
tmp = (sum(MDDM_H + MDDM_V + MDDM_D,3)+MDDM_A)>quantile(vMDDM_boot,.999);
imagesc(tmp)
imagesc(sum(MDDM_H + MDDM_V + MDDM_D,3) + MDDM_A)


vMDDM_boot = randsample(Hel_hat_A,1000,true);
tmp = MDDM_A>quantile(vMDDM_boot,.999);
imagesc(tmp)
imagesc(sum(MDDM_H + MDDM_V + MDDM_D,3))


