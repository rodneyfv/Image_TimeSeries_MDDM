function [ vDimSel, mPvalues, vLamb ] = Estim_Dim_Pval( Xdec, p, Nboot, alpha )
% in this function we test the dimension of curves time series using the
% method of Fonseca and Pinheiro (2019). The function receives a matrix of
% wavelet coefficients of the curves, a lag p and performs bootstrap tests
% with Nboot replications and significance level alpha. The function
% returns the selected dimension vDimSel, the p-values obtained in the
% bootstrap tests and a vector with eigenvalues associanted to the
% eigenfunctions.

[~, n] = size(Xdec);
% vector with mean of the coefficients along time
mu_dec = mean(Xdec,2);

% matrix with the deviations of the coefficients with respect to the mean
% coefficients in a same day
C = Xdec - mu_dec*ones(1,n);

% now we obtain a matrix composed of the coefficients from the decomposition
% of the observed fuctional
C1 = C(:,1:(n-p));
D1 = zeros(n-p,n-p);
for k=1:p
    D1 = D1 + C(:,(k+1):(n-p+k))'*C(:,(k+1):(n-p+k));
end
D = C1*D1*C1'/((n-p)^2);

% after computing D, we obtain its eigenvectors and eigenvalues. The
% columns of B have the eigenvectors of D, which are wavelet coefficients
% of the eigenfunctions corresponding to eigenvalues in the diagonal of L
[B,L] = eig(D);
% reordering the eigenvectors in decreasing order corresponding to the
% eigenvalues
[L,vI] = sort(diag(L),'descend');
B = B(:, vI);
L = diag(L);

vLamb = diag(L);

% performing a bootstrap test to estimate the process' dimension
d0 = 1;
% if the estimated dimension is greater than 8, vDimSel will be
% fixed at 9
vDimSel = 9;
mPvalues = ones(1,8);  % p-values of the tests for the largest 8 eigenvalues

while d0<=8
    % we simulate a process of dimension d0 and then obtain an empirical
    % bootstrap distribution of the (d0+1)-th largest eigenvalue of D using
    % the function DimEst_wavestrap
    d_boot = DimEst_wavestrap( Xdec, Nboot, B(:,[1:d0]), p);
    % p-value for the (d0+1)-th largest eigenvalue of D
    mPvalues(1,d0) = sum(d_boot>L(d0+1,d0+1))/Nboot;
    % we verify now if the null hypothesis that this eigenvalue is zero is
    % rejected or not
    if (mPvalues(1,d0)>alpha)
        % if the null hypothesis is not rejected for d0+1, then the
        % selected dimension is d0
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end



end

