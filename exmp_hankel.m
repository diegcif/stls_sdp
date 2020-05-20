% Compute the nearest rank deficient Hankel matrix

% construct space of m x n hankel matrices
m = 4;
n = 6;
[S,k] = hankel_struct(m,n);

% generate random hankel matrix
u1 = randn(1,k);
u1 = u1/norm(u1);
U1 = applyAffineMap(S,u1);
disp('random hankel matrix')
disp(U1)

% find nearest rank deficient hankel matrix
tic
[opt,u,U,z,X] = sdp_stls(S,u1);
toc
disp('rank deficient hankel matrix')
disp(U)
