function [lamda, vec, error, real, counter, rateOfConvergence] = inverse_shift_power(A)
x = zeros(size(A,1),1); %starting guess, x_0
x(1) = 1;   %to keep x_0 from being equal to zero
x = x/norm(x,2);
sigma = 0;  %no good value found
M = (A - sigma*eye(size(A)));   %shifted matrix M
tolerance = 1e-4;   %required tolerance from assignment
lamda = 0; counter = 0;     %def variables before entering while
lamda_0 = lamda + 2*tolerance;
[L,U,P] = lu(M);    %LUP-composing
while(abs(lamda -lamda_0) > tolerance) 
    counter = counter +1; %no of iterations
    y = U \(L\x);    %iterate until required tolerance is reached
    lamda_0 = lamda;
    lamda = (x'*x)/(x'*y) + sigma;
    [w,i] = max(abs(y));
    alpha = y(i);
    x = y/alpha;    
end
%outputa for function inverse_shift_power
abs(lamda);
real = min(abs(eig(A)));   %real eigenvalue
vec = x;
conv = abs(eig(inv(A)));
error = (abs(lamda-min(eig(A))))/min(abs(eig(A)));
counter;
rateOfConvergence = conv(2)/conv(1);

end
