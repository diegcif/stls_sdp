% Constructs the affine map of Hankel structure
% P: R^n -> R^{k x m}, 
%    u -> hankel(u)

function [PP,n,chordalstr] = hankel_struct(k,m)

n = k+m-1;
Ik = eye(k);

P = zeros(k,n+1,m); % last row of Pi is zero
PP = zeros(k*(n+1),m);
for i=1:m
    P(:,i:i+k-1,i) = Ik;
    PP(:,i) = vect(P(:,:,i));
end

chordalstr = get_struct(k,m,n);

function v = vect(M)
v = M(:);

function chordalstr = get_struct(k,m,n)
M = [hankel(1:k,k:n); (n+1)*ones(1,m)];
cliques = num2cell(M',2)';
eqs = num2cell(1:m);
chordalstr = [cliques; eqs];