% Constructs the affine map of Hankel structure
% P: R^k -> R^{m x n},
%    u -> hankel(u)

function [PP,k] = hankel_struct(m,n)

k = m+n-1;
Im = eye(m);

P = zeros(m,k+1,n); % last row of Pi is zero
PP = zeros(m*(k+1),n);
for i=1:n
    P(:,i:i+m-1,i) = Im;
    PP(:,i) = vect(P(:,:,i));
end

function v = vect(M)
v = M(:);
