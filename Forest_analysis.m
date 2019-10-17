% 

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
J = 5;

% sample size of the time series
n = 87;

% reading the vector with the dates each image was observed. It is saved in
% a cell whose first elements have datetime variables and the second element
% have the image's respective months
tmp1 = readtable('time.txt','ReadVariableNames',false);
cTime = {};
cTime{1} = {};
cTime{2} = zeros(n,1);
for ii=1:n
    tmp2 = char(tmp1{ii,:});
    cTime{1}{ii} = datetime(tmp2(3:12),'InputFormat','yyyy-MM-dd');
    cTime{2}(ii) = month(cTime{1}{ii});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we test the idea of estimating the dimension of the subspace
% generating the density time series of wavelet coefficients, then using
% the densities estimated with the eigenfunctions and the dimension chosen
% to obtain residuals, which are employed to obtain a bootstrap
% distribution of the Hellinger distance under H0.
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we perform the same analysis as Atto el at (2018), but the MDDM
% spots are compared by running a Wilcoxon test on the MDDM's columns to
% identify significantly different time points. Later, we perform the
% wavelet estimation of functional dimension on each density time series,
% and use the loadings to make a one-step ahead prediction of the density
% for each subband on the 2D decomposition.
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

% plotting the sum of the MDDMs with some dates on the x-axis
imagesc(sum(MDDM_H + MDDM_V + MDDM_D,3) + MDDM_A)
tmp = cTime{1};
for ii=1:n
    if(mod(ii,15)==0||ii==2||ii==n)
        tmp{ii} = char(tmp{ii});
    else
        tmp{ii} = '';
    end
end
set(gca,'xtick',[1:n],'xticklabel',tmp)


% matrix with p-values of the Wilcoxon test of the null hypothesis that two
% paired vectors come from a distribution with the same medians. Since
% multiple tests are being made, we set a low critical value for the
% p-values
mMDDMsum = sum(MDDM_H + MDDM_V + MDDM_D,3)+MDDM_A;
mPval = ones(size(mMDDMsum));
for k=1:(n-1)
    for m=(k+1):n
        mPval(k,m) = ranksum(mMDDMsum(:,k),mMDDMsum(:,m));
        mPval(m,k) = mPval(k,m);
    end    
end
imagesc(mPval<1.0e-10)


% arrays to store the mean coefficients of each month
mMonthMean_A = zeros(n_dec,12);
mMonthMean_H = zeros(n_dec,J,12);
mMonthMean_V = zeros(n_dec,J,12);
mMonthMean_D = zeros(n_dec,J,12);
% arrays to store the processes with the monthly means subtracted
Xdec_A_0mean = zeros(n_dec,n);
Xdec_H_0mean = zeros(n_dec,J,n);
Xdec_V_0mean = zeros(n_dec,J,n);
Xdec_D_0mean = zeros(n_dec,J,n);
% computing monthly means and subtracting them from the process'
% coefficients
for t=1:12
    tmp1 = find(cTime{2}==t);
    mMonthMean_A(:,t) = mean(Xdec_A(:,tmp1),2);
    Xdec_A_0mean(:,tmp1) = Xdec_A(:,tmp1) - mMonthMean_A(:,t)*ones(1,length(tmp1));
    for jj=1:J
        mMonthMean_H(:,jj,t) = mean(Xdec_H(:,jj,tmp1),3);
        mMonthMean_V(:,jj,t) = mean(Xdec_V(:,jj,tmp1),3);
        mMonthMean_D(:,jj,t) = mean(Xdec_D(:,jj,tmp1),3);
        Xdec_H_0mean(:,jj,tmp1) = reshape(Xdec_H(:,jj,tmp1),n_dec,length(tmp1)) - mMonthMean_H(:,jj,t)*ones(1,length(tmp1));
        Xdec_V_0mean(:,jj,tmp1) = reshape(Xdec_V(:,jj,tmp1),n_dec,length(tmp1)) - mMonthMean_V(:,jj,t)*ones(1,length(tmp1));
        Xdec_D_0mean(:,jj,tmp1) = reshape(Xdec_D(:,jj,tmp1),n_dec,length(tmp1)) - mMonthMean_D(:,jj,t)*ones(1,length(tmp1));
    end
end

% plot of the monthly densities
% getting the matrix returned by wavedec2
mS = log(SatImages.data(idx,idy,1,1) + 1) - log(img_mean + 1);
[~,~,~,~, L1] = FullDecompImageTS_wavedec2( mS, J, wname, npts, j1 );
% points used for binning of approx. coeff. density
%vpts = linspace(-5,40,npts);
vpts = linspace(-0.5,0.5,npts);
for jj=1:J
    % getting the density estimated at time t in the approximation subband
    %sqrtdens1 = waverec(mMonthMean_H(:,J,jj),L1,wname);
    sqrtdens1 = waverec(mMonthMean_A(:,jj),L1,wname);
    % plotting the densities
    plot(vpts,sqrtdens1.^2,'LineWidth',.5), set(gca,'FontSize',5)
    %ylim([0 16]), title(sprintf('t=%d',t),'FontSize',5)    
    hold on
end

% estimating the dimension of the density time series for each sub-region
% of the 2D decomposition of the images
p = 1;  % maximum lag used
Nboot = 200;  % number of bootstrap replications
alpha = .05;  % significance level of the bootstrap test

dim_A = 0;
vdim_H = ones(J,1);
vdim_V = ones(J,1);
vdim_D = ones(J,1);
[dim_A,~,~] = Estim_Dim_Pval_0mean( reshape(Xdec_A_0mean,[n_dec n]), p, Nboot, alpha );
parfor jj=1:J
    [vdim_H(jj),~,~] = Estim_Dim_Pval_0mean( reshape(Xdec_H_0mean(:,jj,:),[n_dec n]), p, Nboot, alpha );
    [vdim_V(jj),~,~] = Estim_Dim_Pval_0mean( reshape(Xdec_V_0mean(:,jj,:),[n_dec n]), p, Nboot, alpha );
    [vdim_D(jj),~,~] = Estim_Dim_Pval_0mean( reshape(Xdec_D_0mean(:,jj,:),[n_dec n]), p, Nboot, alpha );
end

% obtaining the wavelet coefficients of the eigenfunctions and their
% respective loadings
cEigFunc_H = {}; cLoadings_H = {};
cEigFunc_V = {}; cLoadings_V = {};
cEigFunc_D = {}; cLoadings_D = {};
[ ~, cEigFunc_A, cLoadings_A ] = Estim_Dim_Pval_0mean( reshape(Xdec_A_0mean,[n_dec n]), p, Nboot, alpha, dim_A);
for jj=1:J
    [ ~, cEigFunc_H{jj}, cLoadings_H{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_H_0mean(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_H(jj));
    [ ~, cEigFunc_V{jj}, cLoadings_V{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_V_0mean(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_V(jj));
    [ ~, cEigFunc_D{jj}, cLoadings_D{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_D_0mean(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_D(jj));
end


% Fitting a VAR(3) model for the first two loading vectors
% We regress the loadings on their lagged values
[mBetaEta,mSigmaEta,mResEta,mCovB,~] = mvregress([cLoadings_A(3:(n-1),1:dim_A) cLoadings_A(2:(n-2),1:dim_A) cLoadings_A(1:(n-3),1:dim_A)],...
    cLoadings_A(4:n,1:dim_A));
% Results for mEta1
% The colums have the estimates, standard values, Z-statists and p-values,
% respectively
[mBetaEta(:,1), sqrt(diag(mCovB)),...
    mBetaEta(:,1)./sqrt(diag(mCovB)), 2*normcdf(-abs(mBetaEta(:,1))./sqrt(diag(mCovB)))]
% F-statistic and p-value for the model of mEta1
((n-3-3)/3)*(sum(cLoadings_A(4:n,1).^2) - sum(mResEta(:,1).^2))/sum(mResEta(:,1).^2)
fcdf(((n-3-3)/3)*(sum(cLoadings_A(4:n,1).^2) - sum(mResEta(:,1).^2))/sum(mResEta(:,1).^2), 3, (n-3-3), 'upper')
% one-step-ahead prediction for the loadings
[cLoadings_A(n,1) cLoadings_A((n-1),1) cLoadings_A((n-2),1)]*mBetaEta

% Now we compute one-step ahead prediction of the density for each subband
% of the 2D wavelet decomposition
pred_Loadings_H = {};
pred_Loadings_V = {};
pred_Loadings_D = {};

% Fitting a VAR(3) model and predicting loadings of approximation
[mBetaEta,mSigmaEta,mResEta,mCovB,~] = mvregress([cLoadings_A(3:(n-1),1:dim_A) cLoadings_A(2:(n-2),1:dim_A) cLoadings_A(1:(n-3),1:dim_A)],...
    cLoadings_A(4:n,1:dim_A));
pred_Loadings_A = [cLoadings_A(n,1) cLoadings_A((n-1),1) cLoadings_A((n-2),1)]*mBetaEta;
for jj=1:J
    % horizontal coefficients
    [mBetaEta,mSigmaEta,mResEta,mCovB,~] = mvregress([cLoadings_H{jj}(3:(n-1),1:vdim_H(jj)) cLoadings_H{jj}(2:(n-2),1:vdim_H(jj)) cLoadings_H{jj}(1:(n-3),1:vdim_H(jj))],...
        cLoadings_H{jj}(4:n,1:vdim_H(jj)));
    pred_Loadings_H{jj} = [cLoadings_H{jj}(n,1:vdim_H(jj)) cLoadings_H{jj}((n-1),1:vdim_H(jj)) cLoadings_H{jj}((n-2),1:vdim_H(jj))]*mBetaEta;    
    % vertical coefficients
    [mBetaEta,mSigmaEta,mResEta,mCovB,~] = mvregress([cLoadings_V{jj}(3:(n-1),1:vdim_V(jj)) cLoadings_V{jj}(2:(n-2),1:vdim_V(jj)) cLoadings_V{jj}(1:(n-3),1:vdim_V(jj))],...
        cLoadings_V{jj}(4:n,1:vdim_V(jj)));
    pred_Loadings_V{jj} = [cLoadings_V{jj}(n,1:vdim_V(jj)) cLoadings_V{jj}((n-1),1:vdim_V(jj)) cLoadings_V{jj}((n-2),1:vdim_V(jj))]*mBetaEta;
    % horizontal coefficients
    [mBetaEta,mSigmaEta,mResEta,mCovB,~] = mvregress([cLoadings_D{jj}(3:(n-1),1:vdim_D(jj)) cLoadings_D{jj}(2:(n-2),1:vdim_D(jj)) cLoadings_D{jj}(1:(n-3),1:vdim_D(jj))],...
        cLoadings_D{jj}(4:n,1:vdim_D(jj)));
    pred_Loadings_D{jj} = [cLoadings_D{jj}(n,1:vdim_D(jj)) cLoadings_D{jj}((n-1),1:vdim_D(jj)) cLoadings_D{jj}((n-2),1:vdim_D(jj))]*mBetaEta;
end


% predicted wavelet coefficients related to densities of each subband of
% the 2D decomposition
pred_Xdec_H = zeros(n_dec,J);
pred_Xdec_V = zeros(n_dec,J);
pred_Xdec_D = zeros(n_dec,J);
pred_Xdec_A = cEigFunc_A*pred_Loadings_A + mMonthMean_A(:,12);
for jj=1:J
    pred_Xdec_H(:,jj) = cEigFunc_H{jj}*pred_Loadings_H{jj} + mMonthMean_H(:,jj,12);
    pred_Xdec_V(:,jj) = cEigFunc_V{jj}*pred_Loadings_V{jj} + mMonthMean_V(:,jj,12);
    pred_Xdec_D(:,jj) = cEigFunc_D{jj}*pred_Loadings_D{jj} + mMonthMean_D(:,jj,12);
end

% comparing the Hellinger distance of observed values with the predited
% densities' coefficients
pred_MDDM_A = zeros(n,1);
pred_MDDM_V = zeros(n,J);
pred_MDDM_H = zeros(n,J);
pred_MDDM_D = zeros(n,J);
for k=1:n
    pred_MDDM_A(k) = norm(Xdec_A(:,k) - pred_Xdec_A)/sqrt(2);
    for jj=1:J
        pred_MDDM_H(k,jj) = norm(Xdec_H(:,jj,k) - pred_Xdec_H(:,jj))/sqrt(2);
        pred_MDDM_V(k,jj) = norm(Xdec_V(:,jj,k) - pred_Xdec_V(:,jj))/sqrt(2);
        pred_MDDM_D(k,jj) = norm(Xdec_D(:,jj,k) - pred_Xdec_D(:,jj))/sqrt(2);
    end
end
% plot of the differences
plot(pred_MDDM_A + sum(pred_MDDM_H + pred_MDDM_V + pred_MDDM_D,2))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we run an idea similar to a moving average, where we consider as
% observed images the average of three consecutive images, then we run the
% usual MDDM analysis with Wilcoxon test of the columns corresponding to
% different time points
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

% since we analazy the average of three images, the sample size will
% decrease by 2
n = 85;
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
    im = SatImages.data(idx,idy,1,t:(t+2));
    % Normalizing the log-images before performing the thresholding
    mS = log(sum(im,4)/3 + 1) - log(img_mean + 1);
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

imagesc(sum(MDDM_H + MDDM_V + MDDM_D,3)+MDDM_A)

% matrix with p-values of the Wilcoxon test of the null hypothesis that two
% paired vectors come from a distribution with the same medians. Since
% multiple tests are being made, we set a low critical value for the
% p-values
mMDDMsum = sum(MDDM_H + MDDM_V + MDDM_D,3)+MDDM_A;
mPval = ones(size(mMDDMsum));
for k=1:(n-1)
    for m=(k+1):n
        mPval(k,m) = ranksum(mMDDMsum(:,k),mMDDMsum(:,m));
        mPval(m,k) = mPval(k,m);
    end    
end
imagesc(mPval<1.0e-10)

