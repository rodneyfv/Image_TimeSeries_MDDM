close all
clear all

% Loading the data set
load Lai2005fig4.dat
y = Lai2005fig4(:, 5)'; % y corresponds to the data to be studied
n = length(y); % n corresponds to the sample size
muG = mean(y(y>2)); %average of the elements greater than one
muL = mean(y(y<2)); %average of the elements less than one

% Applying the transformation to take the problem to a regression context
w = (y - muL)/(muG - muL);

% Arguments
s = 1;
delt = 1e-4; % Length of the sub-intervals in the Trapezoidal rule.
min_pts = 2; % min_pts corresponds to the minimum number of sub-interval less one in the Trapezoidal rule.
wJ = ceil(.75*log2(n)); % Considering the resolution level 6
wfilt = [0.027333068345164, 0.029519490926073,-0.039134249302583,...
         0.199397533976996, 0.723407690403808, 0.633978963456949,...
         0.016602105764424,-0.175328089908107,-0.021101834024930,...
         0.019538882735387]; %Daubechies Least Asymmetric 10-tap wavelet filter called wfilter
wprec = 30;

% We consider as the time of observation a grid of points in the unit
% interval
x = (1:n)/n; 

% Obtaining the raw scaling coefficient estimates
[w1, w2] = WavCoefEsts(w, x, s, delt, min_pts, wJ, wfilt, wprec);

% Decomposing the raw estimates into wavelet domain
[wd1, snoise1] = wd(w1, 0, wfilt, true);
wdH1 = wd1(1:2^wJ);
[wd2, snoise2] = wd(w2, 0, wfilt, true);
wdH2 = wd2(1:2^wJ);

% Obtaining the raw function estimates
yy = PHI(x, wJ, wfilt, wprec, true);
f1row1 = yy*w1';
f2row1 = yy*w2';

% Shrinkage using hard threshold
rwdH1 = WavShrink(wdH1, wfilt, true, 'Hard', true, snoise1);
rwdH2 = WavShrink(wdH2, wfilt, true, 'Hard', true, snoise2);

% Function estimates regularized by Hard threshold
fH1 = yy*rwdH1';
fH1(fH1<0) = 0;
fH1(fH1>1) = 1;
fH2 = yy*rwdH2';
fH2(fH2<0) = 0;
fH2(fH2>1) = 1;

% Shrinkage using Efromovich's approach
[rwdEfr1, Jest1, risk1] = WavShrink(wd1, wfilt, true, 'Efromovich', true, snoise1);
[rwdEfr2, Jest2, risk2] = WavShrink(wd2, wfilt, true, 'Efromovich', true, snoise2);

% Function estimates regularized by Efromovich's shrinkage
yy1 = PHI(x, Jest1, wfilt, wprec, true);
yy2 = PHI(x, Jest2, wfilt, wprec, true);
fEfr1 = yy1*rwdEfr1';
fEfr1(fEfr1<0) = 0;
fEfr1(fEfr1>1) = 1;
fEfr2 = yy2*rwdEfr2';
fEfr2(fEfr2<0) = 0;
fEfr2(fEfr2>1) = 1;

% Function estimates regularized by Local Linear Regression
r1=ksrlin(x,f1row1);
r1.f(r1.f<0) = 0;
r1.f(r1.f>1) = 1;
r2=ksrlin(x,f2row1);
r2.f(r2.f<0) = 0;
r2.f(r2.f>1) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ploting the function estimates as in Figure 11
% filename = '~/fEstimates.pdf';

lw = 3;  
set(0, 'DefaultAxesFontSize', 13);
fs = 16;
msize = 10;

h = figure;

subplot(3, 2, 1)
plot(x, fH1', 'k-')
% xlabel('') % hard threshold

subplot(3, 2, 2)
plot(x, fH2', 'k-')
% xlabel('') % hard threshold

subplot(3, 2, 3)
plot(x, fEfr1', 'k-')
% xlabel('') % hard threshold

subplot(3, 2, 4)
plot(x, fEfr2', 'k-')
% xlabel('') % hard threshold

subplot(3, 2, 5)
plot(linspace(0, 1, 250), r1.f, 'k-')
xlabel('Index') % 

subplot(3, 2, 6)
plot(linspace(0, 1, 250), r2.f, 'k-')
xlabel('Index') % hard threshold

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,filename,'-dpdf','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the data set (Figure 10)
% filename = '~/GBMdata.pdf';

% lw = 3;  
% set(0, 'DefaultAxesFontSize', 13);
% fs = 16;
% msize = 10;

h = figure;

scatter(x, y, 'k')
xlabel('Index')

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,filename,'-dpdf','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
