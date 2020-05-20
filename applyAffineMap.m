% Applies the affine map
% P: R^k -> R^{m x n}
% to a vector u in R^k

function U = applyAffineMap(PP,u)

k = length(u);
n = size(PP,2);
m = size(PP,1)/(k+1);

v = [u(:);1];
U = zeros(m,n);
for i=1:n
    U(:,i) = reshape(PP(:,i),[m,k+1])*v;
end
