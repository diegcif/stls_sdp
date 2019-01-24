% Applies the affine map
% P: R^n -> R^{k x m}
% to a vector u in R^n

function U = applyAffineMap(PP,u)

n = length(u);
m = size(PP,2);
k = size(PP,1)/(n+1);

v = [u(:);1];
U = zeros(k,m);
for i=1:m
    U(:,i) = reshape(PP(:,i),[k,n+1])*v;
end