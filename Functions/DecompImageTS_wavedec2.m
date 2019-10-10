function [ xdec_H, xdec_V, xdec_D, L1 ] = DecompImageTS_wavedec2( img_obs, J, wname, npts, j1 )
% 


% changing the current boundary to periodic.
dwtmode('per','nodisp');

% smoothing the image using 2D wavelet denoising
[thr,~,~] = ddencmp('cmp','wp',img_obs);
mG = wdencmp('gbl',img_obs,'db6',6,thr*2,'s',1);

% taking the 2D-DWT of the thresholded normalized image
[vC,mL] = wavedec2(mG,J,wname);

% getting the matrices with detail coefficients up to level J_MDDM and
% storing them as a vector in a cell
cellH = cell(J);
cellV = cell(J);
cellD = cell(J);
for jj=1:J
    [tmpH,tmpV,tmpD] = detcoef2('all',vC,mL,jj);
    cellH{jj} = tmpH(:);
    cellV{jj} = tmpV(:);
    cellD{jj} = tmpD(:);
end

% using the same points of binning for all density estimates
vpts = linspace(-0.5,0.5,npts);
% window length used in the binning
wind = (vpts(npts) - vpts(1))/sqrt(npts);

[~,C1,~] = wden(ones(npts,1),'sqtwolog','s','sln',j1,wname);
n_dec = length(C1);

xdec_H = zeros(n_dec,J);
xdec_V = zeros(n_dec,J);
xdec_D = zeros(n_dec,J);

for jj=1:J
    % vectorizing the horizontal detail coefficients stored in cellH
    % binning, which gives a pre-estimate of the square-root of the densities
    bin1 = sqrt(binning(cellH{jj},vpts,wind));
    % wavelet decomposition and denoising of the binned data
    [~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    xdec_H(:,jj) = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

    % performing the same for vertical detail coefficients
    bin1 = sqrt(binning(cellV{jj},vpts,wind));
    [~,C1,~] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    xdec_V(:,jj) = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));

    % performing the same for diagonal detail coefficients
    bin1 = sqrt(binning(cellD{jj},vpts,wind));
    [~,C1,L1] = wden(bin1,'sqtwolog','s','sln',j1,wname);
    xdec_D(:,jj) = C1/sqrt(sum(C1.^2)*(vpts(2)-vpts(1)));
end


end

