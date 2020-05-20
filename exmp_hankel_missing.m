% Nearest rank deficient Hankel matrix with missing entries

% construct space of m x n hankel matrices
m = 4;
n = 6;
[S,k] = hankel_struct(m,n);

% missing indices
i_miss = [m+1,k-m+1];

% generate random hankel matrix
u1 = randn(1,k);
u1(i_miss) = nan;
U1 = hankel(u1(1:m),u1(m:k));
disp('random hankel matrix')
disp(U1)

% find nearest rank deficient hankel matrix
tic
[opt,u,U,z,X] = sdp_stls(S,u1);
toc
disp('rank deficient hankel matrix')
disp(U)
