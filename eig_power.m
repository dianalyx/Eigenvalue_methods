function [eigenvalue,realvalue,error,iterations, rateOfConvergence, eigenvector] = eig_power(A)
x = ones(size(A,1),1);
tolerance = 1e-4;
x = x/norm(x,2);
lamda = 0;
lamda_0 = lamda + 2*tolerance;
counter = 0;
y = [];
while(abs(lamda -lamda_0) > tolerance)
    y = A*x;
    lamda_0 = lamda;
    lamda = x'*y;
    x = y/norm(y,2);
    counter = counter+1;
end
eigenvalue = lamda;
eigenvector = y;
realvalue = abs(max(eig(A)));
conv = abs(eig(A));
error = realvalue - eigenvalue; %Calculates the
error = error/realvalue;        %error in percentage.
iterations = counter;
rateOfConvergence = conv(2)/conv(1);
plot(1:length(A),eigenvector);
xlabel('Eigenvector');
ylabel('Value of normalized eigenvector');
end