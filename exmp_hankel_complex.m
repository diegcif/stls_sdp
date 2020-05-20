% Nearest rank deficient Hankel matrix with complex data

% construct space of m x n hankel matrices
m = 3;
n = 3;
[S,k] = hankel_struct(m,n);

% generate random hankel matrix
u1 = randn(k,1) + 1i*randn(k,1);
U1 = applyAffineMap(S,u1);
disp('random hankel matrix')
disp(U1)

% find nearest rank deficient hankel matrix
tic
[opt,u,U,z,X] = sdp_stls(S,u1);
toc
disp('rank deficient hankel matrix')
disp(U)