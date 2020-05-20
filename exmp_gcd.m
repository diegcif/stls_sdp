% Solve approximate GCD problem

% construct space of sylvester matrices
D = [6 5]; % degrees of f1 and f2
d = 4; % degree of the gcd
noise = 5e-1;
S = sylvester_struct(D,d);

% generate random polynomials
g = randn(d+1,1); % gcd
h1 = randn(D(1)-d+1,1);
h2 = randn(D(2)-d+1,1);
f1 = conv(g,h1); f1 = f1/norm(f1); % poly1
f2 = conv(g,h2); f2 = f2/norm(f2); % poly2
f1 = f1 + noise*randn(D(1)+1,1); % add noise
f2 = f2 + noise*randn(D(2)+1,1);

% construct Sylvester matrix
u1 = [f1; f2];
U1 = applyAffineMap(S,u1);
disp('random Sylvester matrix')
disp(U1)

% find nearest rank deficient Sylvester matrix
tic
[opt,u,U,z,X] = sdp_stls(S,u1);
toc
disp(opt)
disp('rank deficient Sylvester matrix')
disp(U)
