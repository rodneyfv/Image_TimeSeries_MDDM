function [ Yort ] = MatGramSchm(Y)
% function to apply the Gram-Schmidt to the vectors of the matrix Y

[n,d] = size(Y);
% matrix to keep the orthonormalized vectors
Yort = zeros(n,d);
Yort(:,1) = Y(:,1);

% steps of orthogonalization
for ii=2:d
    Yort(:,ii) = Y(:,ii);
    for jj=1:(ii-1)
        Yort(:,ii) = Yort(:,ii) - ((Yort(:,jj)'*Y(:,ii))/(Yort(:,jj)'*Yort(:,jj)))*Yort(:,jj);
    end
end

% steps to normalization
for ii=1:d
    Yort(:,ii) = Yort(:,ii)/norm(Yort(:,ii));
end

end

