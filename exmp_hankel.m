% Compute the nearest rank deficient Hankel matrix

% construct space of k x m hankel matrices
k = 4;
m = 6;
[PP,n] = hankel_struct(k,m);

% generate random hankel matrix
u1 = randn(1,n);
u1 = u1/norm(u1);
U1 = applyAffineMap(PP,u1);
disp('random hankel matrix')
disp(U1)

% find nearest rank deficient hankel matrix
tic
[opt,u,U,z,X] = sdp_stls(PP,u1);
toc
disp('rank deficient hankel matrix')
disp(U)