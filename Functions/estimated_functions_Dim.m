function [ mY ] = estimated_functions_Dim( Xdec, p, d )
% this function receives a matrix Xdec with wavelet coefficients and using
% a lag of p, it returns the coefficients mY of a estimated process of
% dimension d.

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
[~,vI] = sort(diag(L),'descend');
B = B(:, vI);

% obtaining wavelet coefficients corresponding to functions generated from
% a subspace with dimension d
mY = mu_dec*ones(1,n);
for ii=1:d
    mY = mY + B(:,ii)*(C'*B(:,ii))';
end


end

