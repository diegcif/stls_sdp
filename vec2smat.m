function M = vec2smat(v)

N = length(v);
n = (-1+sqrt(1+8*N))/2;

M = repmat(0*v(1:n),[1,n]);
I = triu(true(n,n),0);
I2 = triu(true(n,n),1);
M(I) = v/2;
M(I2) = M(I2)*sqrt(2);
M = M + M.';
