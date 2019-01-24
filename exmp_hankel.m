% hankel matrices of size k x m
k = 2;
m = 3;
[PP,n,chordalstr] = hankel_struct(k,m);

% random vector
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