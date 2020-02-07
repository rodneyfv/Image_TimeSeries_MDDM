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

% making a loop to generate txt files that will be used in R to run a
% spatial model in the images
for ii=1:87
    img = log(SatImages.data(idx,idy,1,ii) + 1);
    dlmwrite(sprintf( 'ImagesTStxt/image-%i.txt',ii),img);
end

% and this is how we can read a treated image
img2 = dlmread('imagehs-1.txt');


