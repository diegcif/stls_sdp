% Constructs the affine map of Sylvester structure
% P: R^D -> R^{k x m},
%    (u1,u2) -> Syl(u1,u2)

function [PP,n,k,m] = sylvester_struct(D,d)
n = sum(D+1);
k = sum(D)-2*d+2;
m = sum(D)-d+1;

PP = zeros(k*(n+1),m);
for i=1:sum(D)-d+1
    Pi = get_matrixGCD_i(D,d,i);
    Pi = [Pi zeros(k,1)];
    PP(:,i) = Pi(:);
end

function Pi = get_matrixGCD_i(D,d,i)
k1 = D(2)-d+1;
k2 = D(1)-d+1;
P1 = (eye(D(1)+1, D(1)+1));
P1 = [zeros(D(2)-d,D(1)+1); P1; zeros(D(2)-d,D(1)+1)];
P2 = P1(end-i+2-k1:end-i+1,:);
Q1 = (eye(D(2)+1, D(2)+1));
Q1 = [zeros(D(1)-d,D(2)+1); Q1; zeros(D(1)-d,D(2)+1)];
Q2 = Q1(end-i+2-k2:end-i+1,:);
Pi = blkdiag(P2,Q2);