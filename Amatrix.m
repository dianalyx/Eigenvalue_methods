clear;
%construct matrix
%n = 20,40,80,200
n = 40;
h = 1/n;

K1 = 2*diag(diag(ones(n)));
K2 = -1*diag(diag(ones(n-1)),1);
K3 = -1*diag(diag(ones(n-1)),-1);
K = K3+K1+K2;
K(end,end) = 1;
K = (1/h)*K;

M1 = 4*diag(diag(ones(n)));
M2 = 1*diag(diag(ones(n-1)),-1);
M3 = 1*diag(diag(ones(n-1)),1);
M = M3+M1+M2;
M(end,end) = 2;
M = (h/6)*M;

A = inv(M)*K;