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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this part we read the images, perform the log transform and then save
% them in text files to read later in R, where the kriking will be
% performed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % making a loop to generate txt files that will be used in R to run a
% % spatial model in the images
% for ii=1:87
%     img = log(SatImages.data(idx,idy,1,ii) + 1);
%     dlmwrite(sprintf( 'ImagesTStxt/image-%i.txt',ii),img);
% end
% 
% % and this is how we can read a treated image
% img2 = dlmread('imagehs-1.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the ACFs of the pixels time series from images smoothed with
% kriging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imagesTS = zeros(512,512,n);
for t=1:87
    % reading the (log) image smoothed with kriging
    imagesTS(:,:,t) = dlmread(fullfile(pwd(),'..','Images_TimeSeries_files/ImagesHatTS',sprintf('imagehs-%i.txt',t)));
end

% estimating the autocorrelation function of the pixel time series
acfsTS = zeros(512,512,21);
vMeanACF = zeros(21,1);  % mean of autocorrelation functions
for ii=1:512
    for jj=1:512
        acfsTS(ii,jj,:) = autocorr(imagesTS(ii,jj,:));
        vMeanACF = vMeanACF + squeeze(acfsTS(ii,jj,:));
    end
    if mod(round(100*ii/512),20)==0;fprintf('%d%%,',round(100*ii/512));end
end
vMeanACF = vMeanACF./(512*512);

% ploting some autocorrelation functions together with the mean in a
% different color
plot(squeeze(acfsTS(1,1,:)),'b')
ylim([-.6 1])
hold on
tic
for ii=round(linspace(1,512,50))
    for jj=round(linspace(1,512,50))
        plot(squeeze(acfsTS(ii,jj,:)),'b')
    end
end
toc
plot(vMeanACF,'r','LineWidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the ACFs of the pixels time series from the difference of images
% smoothed with kriking and the observed images, that is, of noise images
% time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoiseTS = zeros(512,512,n);
for t=1:87
    % image smoothed with kriking
    imgsm = dlmread(fullfile(pwd(),'..','Images_TimeSeries_files/ImagesHatTS',sprintf('imagehs-%i.txt',t)));
    % raw (log) image
    imgraw = dlmread(fullfile(pwd(),'..','Images_TimeSeries_files/ImagesTStxt',sprintf('image-%i.txt',t)));
    % noise from this smoothing
    NoiseTS(:,:,t) = imgraw - imgsm;
end

% estimating the autocorrelation function of the pixel time series
acfsTS = zeros(512,512,21);
vMeanACF = zeros(21,1);  % mean of autocorrelation functions
for ii=1:512
    for jj=1:512
        acfsTS(ii,jj,:) = autocorr(NoiseTS(ii,jj,:));
        vMeanACF = vMeanACF + squeeze(acfsTS(ii,jj,:));
    end
    if mod(round(100*ii/512),20)==0;fprintf('%d%%,',round(100*ii/512));end
end
vMeanACF = vMeanACF./(512*512);

% ploting some autocorrelation functions together with the mean in a
% different color
plot(squeeze(acfsTS(1,1,:)),'b')
ylim([-.6 1])
hold on
tic
for ii=round(linspace(1,512,50))
    for jj=round(linspace(1,512,50))
        plot(squeeze(acfsTS(ii,jj,:)),'b')
    end
end
toc
plot(vMeanACF,'r','LineWidth',2)


