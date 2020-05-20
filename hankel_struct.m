% Constructs the affine map of Hankel structure
% SS: R^k -> R^{m x n},
%    u -> hankel(u)

function [S,k] = hankel_struct(m,n)

k = m+n-1;
Im = eye(m);

SS = zeros(m,k+1,n);
S = zeros(m*(k+1),n);
for i=1:n
    SS(:,i:i+m-1,i) = Im;
    S(:,i) = vect(SS(:,:,i));
end

function v = vect(M)
v = M(:);
