% 

close all
clear all

addpath(genpath('Functions/'))
addpath(genpath('SJS_WBEMR/'))

SatImages = matfile('/media/rodney/Arquivos/Datasets/Abdou_87_Sattelite_Images/xSpatial_ySpatial_VH_VV_Time.mat');
% SatImages.data is a matrix of dimension 4000x3000x2x87

ds = 2;
idx = 1:ds:1024;
idy = 501:ds:(500 + 1024);

% img = SatImages.data(idx,idy,1,1);
% dlmwrite('image0-0.txt',img)

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
% number of coefficients in the wavelet decomposition of the functions
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

for jj=1:J
    % estimated density
    Xdec_hat_H(:,jj,:) = estimated_functions_Dim( reshape(Xdec_H(:,jj,:),[npts n]), p, vdim_H(jj) );
    Xdec_hat_V(:,jj,:) = estimated_functions_Dim( reshape(Xdec_V(:,jj,:),[npts n]), p, vdim_V(jj) );
    Xdec_hat_D(:,jj,:) = estimated_functions_Dim( reshape(Xdec_D(:,jj,:),[npts n]), p, vdim_D(jj) );
end
Xdec_hat_A = estimated_functions_Dim( reshape(Xdec_A,[npts n]), p, dim_A );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we apply the method of Montoril et al (2019) for the unidimensional
% time series of loadings estimated above. For each time series, a mixture
% function is estimated, were we can observe which points seem to be
% significantly different
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtaining the wavelet coefficients of the eigenfunctions and their
% respective loadings
cEigFunc_H = {}; cLoadings_H = {};
cEigFunc_V = {}; cLoadings_V = {};
cEigFunc_D = {}; cLoadings_D = {};
[ ~, cEigFunc_A, cLoadings_A ] = Estim_Dim_Pval_0mean( reshape(Xdec_hat_A,[n_dec n]), p, Nboot, alpha, dim_A);
for jj=1:J
    [ ~, cEigFunc_H{jj}, cLoadings_H{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_hat_H(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_H(jj));
    [ ~, cEigFunc_V{jj}, cLoadings_V{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_hat_V(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_V(jj));
    [ ~, cEigFunc_D{jj}, cLoadings_D{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_hat_D(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_D(jj));
end

% Arguments
s = 1; % Value used to define the intervals
delt = 1e-4; % Length of the sub-intervals in the Trapezoidal rule.
min_pts = 2; % min_pts corresponds to the minimum number of sub-interval less one in the Trapezoidal rule.
wJ = ceil(.75*log2(n)); % Considering the resolution level 6
wfilt = [0.027333068345164, 0.029519490926073,-0.039134249302583,...
         0.199397533976996, 0.723407690403808, 0.633978963456949,...
         0.016602105764424,-0.175328089908107,-0.021101834024930,...
         0.019538882735387]; %Daubechies Least Asymmetric 10-tap wavelet filter called wfilter
wprec = 30; % number of Daubechies-Lagarias steps

% We consider as the time of observation a grid of points in the unit
% interval
x = (1:n)/n; 
rawest = 'wavcoef';
estimator = 'Efromovich';

% estimating the mixture function for approximation loadings
tmp = zeros(sum(dim_A), n);
cont = 1;
for ii=1:dim_A
    tmp(cont,:) = mixtureProbs(cLoadings_A(:,ii), x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator);
    cont = cont + 1;
end
% mean mixture function estimated
vMixtureProbs_A = mean(tmp,1);

% estimating the mixture function for detail loadings separately for each
% detail level
mMixtureProbs_details = zeros(J,n);
for jj=1:J
    tmp = zeros(sum(vdim_H(jj) + vdim_V(jj) + vdim_D(jj)), n);
    cont = 1;
    for ii=1:vdim_H(jj)
        tmp(cont,:) = mixtureProbs(cLoadings_H{jj}(:,ii), x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator);
        cont = cont + 1;
    end
    for ii=1:vdim_V(jj)
        tmp(cont,:) = mixtureProbs(cLoadings_V{jj}(:,ii), x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator);
        cont = cont + 1;
    end
    for ii=1:vdim_D(jj)
        tmp(cont,:) = mixtureProbs(cLoadings_D{jj}(:,ii), x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator);
    end
    mMixtureProbs_details(jj,:) = mean(tmp,1);
end

% plot of the mean mixture function of loadings corresponding to
% approximation and detail coefficients
vMeanMixtureFunc = mean([vMixtureProbs_A; mMixtureProbs_details(1,:)],1);
plot(1:n,vMeanMixtureFunc)
xlabel('t')
ylabel('\rho(t)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we apply the 2D DWT for the images and later apply the method of
% Montoril et al (2019) on the coefficients time series. The mixture
% function obtained for each coeffient is compared with the mixture
% function obtained for the whole image time series while computing the
% MDDM. Then, we check wich coeffients have the mixture function closer to
% its corresponding function for the all images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[vC,mL] = wavedec2(img_mean,J,wname);
n_C = length(vC);
% matrix to keep the 2D DWT of the images
mDWT2D_TS = zeros(n,n_C);
parfor jj=1:2; end
parfor t=1:n
    %SatImages = matfile('/media/rodney/Arquivos/Datasets/Abdou_87_Sattelite_Images/xSpatial_ySpatial_VH_VV_Time.mat');
    % changing the current boundary to periodic.
    dwtmode('per','nodisp');
    % Normalizing the log-images before performing the thresholding
    mS = log(SatImages.data(idx,idy,1,t) + 1) - log(img_mean + 1);
    
    % smoothing the image using 2D wavelet denoising
    [thr,~,~] = ddencmp('cmp','wp',mS);
    mG = wdencmp('gbl',mS,'db6',6,thr*2,'s',1);
    % taking the 2D-DWT of the thresholded normalized image
    [mDWT2D_TS(t,:),~] = wavedec2(mG,J,wname);
end

% Indexes of the coefficients whose time series we'll analyze.
% prod(mL(1,:)) is the number of approximation coefficients in the 2D-DWT,
% and the coefficients we take below correspond to them
vC_ids = 1:prod(mL(1,:));
% matrix that keeps the mixture function of each wavelet coefficient
mMixtureFuncCoef = zeros(n,length(vC_ids));
% vector to keep the norm of the difference of each coefficient mixture
% function with the mean mixture function computed for the whole time
% series
vNormDiff_C = zeros(length(vC_ids),1);
tic
for ii=1:length(vC_ids)
   mMixtureFuncCoef(:,ii) = mixtureProbs(mDWT2D_TS(:,vC_ids(ii)), x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator); 
   vNormDiff_C(ii) = norm(mMixtureFuncCoef(:,ii) - vMeanMixtureFunc');
   %fprintf('ii=%i\n',ii)
end
toc
% index corresponding to the coefficient with lowest norm found
[tmp1, tmp2] = min(vNormDiff_C)
vC_ids(tmp2)
% comparing the mean mixture function with this coefficient mixture
% function
plot(mMixtureFuncCoef(:,tmp2))      
hold on
plot(vMeanMixtureFunc)

% identifying the pixels related to the 10 coefficients with lowest values
% of norm difference
[tmp1, tmp2] = sort(vNormDiff_C,'ascend');
vC = zeros(n_C,1);
vC(vC_ids(tmp2(1:10))) = 1;
tmp = waverec2(vC,mL,wname);
% plotting these points above the mean image
imshow(img_mean)
hold on
[tmp1, tmp2] = find(abs(tmp)>.01);
plot(tmp2, tmp1, 'r*', 'LineWidth', .1, 'MarkerSize', .1);

% Making a video where we plot the change regions chaging while we change
% the proportion of used lowest norms corresponding to approximation
% coefficients.
% create the video writer with 1 fps
writer_im = VideoWriter('ChangingRegions.avi');
% number of frames shown in a second
writer_im.FrameRate = 1;
% open the video writer
open(writer_im);
% write the frames to the video
imshow(img_mean) % plotting the mean image for reference
for t = round(linspace(1,256,20))
    [tmp1, tmp2] = sort(vNormDiff_C,'ascend');
    vC = zeros(n_C,1);
    vC(vC_ids(tmp2(1:t))) = 1;
    tmp = waverec2(vC,mL,wname);
    imshow(img_mean) % plotting the mean image for reference
    hold on
    % plotting these points above the mean image
    [tmp1, tmp2] = find(abs(tmp)>.01);
    plot(tmp2, tmp1, 'r*', 'LineWidth', .1, 'MarkerSize', .1);
    title(sprintf('%d%%',round(100*t/256))) 
    hold off
    frame = getframe(gcf);
    writeVideo(writer_im, frame);    
end
% close the writer object
close(writer_im);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% making a video with the image time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create the video writer with 1 fps
writer_im = VideoWriter('ForestTS.avi');
% number of frames shown in a second
writer_im.FrameRate = 6;
% open the video writer
open(writer_im);
% write the frames to the video
for t = 2:87
    % changing the current boundary to periodic.
    dwtmode('per','nodisp');
    % Normalizing the log-images before performing the thresholding
    mS = log(SatImages.data(idx,idy,1,t) + 1);
    
    % smoothing the image using 2D wavelet denoising
    [thr,~,~] = ddencmp('cmp','wp',mS);
    mG = wdencmp('gbl',mS,'db6',6,thr*2,'s',1);
    
    imshow(exp(mG) - 1)
    title(sprintf('%c',char(cTime{1}{t})))
    frame = getframe(gcf);
    writeVideo(writer_im, frame);    
end
% close the writer object
close(writer_im);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we apply the method of Montoril et al (2019) for the unidimensional
% time series of loadings estimated above. Before applying it, we use a
% k-means algorithm to separate the observations in two groups, which are
% used to compute the different means necessary to apply Montoril's idea.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtaining the wavelet coefficients of the eigenfunctions and their
% respective loadings
cEigFunc_H = {}; cLoadings_H = {};
cEigFunc_V = {}; cLoadings_V = {};
cEigFunc_D = {}; cLoadings_D = {};
[ ~, cEigFunc_A, cLoadings_A ] = Estim_Dim_Pval_0mean( reshape(Xdec_hat_A,[n_dec n]), p, Nboot, alpha, dim_A);
for jj=1:J
    [ ~, cEigFunc_H{jj}, cLoadings_H{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_hat_H(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_H(jj));
    [ ~, cEigFunc_V{jj}, cLoadings_V{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_hat_V(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_V(jj));
    [ ~, cEigFunc_D{jj}, cLoadings_D{jj} ] = Estim_Dim_Pval_0mean( reshape(Xdec_hat_D(:,jj,:),[n_dec n]), p, Nboot, alpha, vdim_D(jj));
end

% stacking the loadings in a single matrix
mLoadings = zeros(dim_A + sum(vdim_H + vdim_V + vdim_D),n);
tmp = dim_A;
mLoadings(1:tmp,:) = cLoadings_A';
tmp = tmp + 1;
for jj=1:J
    for ii=1:vdim_H(jj)
        mLoadings(tmp,:) = cLoadings_H{jj}(:,ii);
        tmp = tmp + 1;
    end
    for ii=1:vdim_V(jj)
        mLoadings(tmp,:) = cLoadings_V{jj}(:,ii);
        tmp = tmp + 1;
    end
    for ii=1:vdim_D(jj)
        mLoadings(tmp,:) = cLoadings_D{jj}(:,ii);
        tmp = tmp + 1;
    end
end
% applying the k-means to separate each of the n observations in two groups
vMixtIdx = kmeans(mLoadings',2);
vMixtIdx = vMixtIdx - 1;

% Arguments
s = 1; % Value used to define the intervals
delt = 1e-4; % Length of the sub-intervals in the Trapezoidal rule.
min_pts = 2; % min_pts corresponds to the minimum number of sub-interval less one in the Trapezoidal rule.
wJ = ceil(.75*log2(n)); % Considering the resolution level 6
wfilt = [0.027333068345164, 0.029519490926073,-0.039134249302583,...
         0.199397533976996, 0.723407690403808, 0.633978963456949,...
         0.016602105764424,-0.175328089908107,-0.021101834024930,...
         0.019538882735387]; %Daubechies Least Asymmetric 10-tap wavelet filter called wfilter
wprec = 30; % number of Daubechies-Lagarias steps

% We consider as the time of observation a grid of points in the unit
% interval
x = (1:n)/n; 
rawest = 'wavcoef';
estimator = 'Efromovich';

% estimating the mixture function for approximation loadings
tmp = zeros(sum(dim_A), n);
cont = 1;
for ii=1:dim_A
    tmp(cont,:) = mixtureProbs_group(normalize(cLoadings_A(:,ii)),vMixtIdx, x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator);
    cont = cont + 1;
end
% mean mixture function estimated
vMixtureProbs_A = mean(tmp,1);

% estimating the mixture function for detail loadings separately for each
% detail level
mMixtureProbs_details = zeros(J,n);
for jj=1:J
    tmp = zeros(sum(vdim_H(jj) + vdim_V(jj) + vdim_D(jj)), n);
    cont = 1;
    for ii=1:vdim_H(jj)
        tmp(cont,:) = mixtureProbs_group(normalize(cLoadings_H{jj}(:,ii)),vMixtIdx, x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator);
        cont = cont + 1;
    end
    for ii=1:vdim_V(jj)
        tmp(cont,:) = mixtureProbs_group(normalize(cLoadings_V{jj}(:,ii)),vMixtIdx, x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator);
        cont = cont + 1;
    end
    for ii=1:vdim_D(jj)
        tmp(cont,:) = mixtureProbs_group(normalize(cLoadings_D{jj}(:,ii)),vMixtIdx, x, s, delt, min_pts, wJ, wfilt, wprec,rawest,estimator);
    end
    mMixtureProbs_details(jj,:) = mean(tmp,1);
end

% plot of the mean mixture function of loadings corresponding to
% approximation and detail coefficients
vMeanMixtureFunc = mean([vMixtureProbs_A; mMixtureProbs_details(1,:)],1);
plot(1:n,vMeanMixtureFunc)
xlabel('t')
ylabel('\rho(t)')

