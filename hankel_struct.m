% Constructs the affine map of Hankel structure
% P: R^n -> R^{k x m}, 
%    u -> hankel(u)

function [PP,n] = hankel_struct(k,m)

n = k+m-1;
Ik = eye(k);

P = zeros(k,n+1,m); % last row of Pi is zero
PP = zeros(k*(n+1),m);
for i=1:m
    P(:,i:i+k-1,i) = Ik;
    PP(:,i) = vect(P(:,:,i));
end

function v = vect(M)
v = M(:);