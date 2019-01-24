% sylvester structure
D = [6 6]; %degf1 degf2
d = 4; % deg gcd
[PP,n,k,m] = sylvester_struct(D,d);

% random vector
f1 = randn(1,D(1)+1); f1 = f1/norm(f1);
f2 = randn(1,D(2)+1); f2 = f2/norm(f2);

f1 = conv([1 -3 2],[1 1]) + [0 0 0 .01];
f2 = conv([1 -3 2],[1 1.2]) - [0 0 0 .01];

f1 = conv(conv([-3 2],[2 1]),[8 -3]);
f2 = conv(conv([-3.02 1.99],[1.998 .98]),[4.5 1]);

g = randn(1,d+1);
% g = [1 2 1];
h1 = randn(1,D(1)-d+1);
h2 = randn(1,D(2)-d+1);
% h1 = [1 3 -2 1];
% h2 = [3 -1 1 -2];
f1 = conv(g,h1); f1 = f1/norm(f1);
f2 = conv(g,h2); f2 = f2/norm(f2);
f1 = f1 + 5e-1*randn(1,D(1)+1);
f2 = f2 + 5e-1*randn(1,D(2)+1);

u1 = [f1 f2];
U1 = applyAffineMap(PP,u1);
% disp('random Sylvester matrix')
% disp(U1)

% find nearest rank deficient Sylvester matrix
tic
[opt,u,U,z,X] = sdp_stls(PP,u1);
toc
disp(opt)
% disp('rank deficient Sylvester matrix')
% disp(U)