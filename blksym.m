% Block symmetrization of a matrix
% Input: mk x mk matrix A
% Output: symmetric matrix S such that 
%         each k x k block is also symmetric (there are m^2 blocks)

function S = blksym(m,k,A)

S = mysym(A);
if isempty(A); return; end

for i=1:m
    for j=1:m
        I = k*(i-1) + (1:k);
        J = k*(j-1) + (1:k);
        S(I,J,:) = mysym(S(I,J,:));
    end
end

function As = mysym(A)
At = permute(A,[[2,1],3:ndims(A)]);
As = .5*(A+At);