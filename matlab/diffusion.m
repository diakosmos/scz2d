clear;
n=100;
g=[0,1,0;0,1,0];
[x,d2,d1,xm,xp]=cheb2bc(n,g);
%[x,d] = chebdif(n,2);
%d1 = d(:,:,1); d2 = d(:,:,2);
a = 0.4;  % +- zero of derivative of function below - other zeros of derivative are +-1.
f= 3*x.^5 - 5*(a^2 + 1)*x.^3 + 15 * a^2 .* x;  % polynomial that is symmetric, nontrivial, and
       % satisfies boundary conditions
f = f + (1 - x.^2);   % make its integral nonzero, to test conservation, since integral of f should
% be conserved, with bc that f's derivative is zero on boundaries.
base_sum_1 = sum(f);       
base_sum_2 = sum(f.^2);
hold off;
plot(x,f);       
    hold on;   
iterations = 1000; i=0;
dt = 0.01;
while i < iterations
    i = i + 1;
    base_sum(i) = sum(f);
    f = f ./ (base_sum(i) / base_sum_1);
    plot(x,f);
    
    %pause;
    
    %d2f = d2 * f;

    %dtf = d2f * dt;

    %f = f + dtf;
    
    mymat = eye(n) - dt * d2;
    %invmymat = 1/mymat;
    f = mymat \ f;
   % temp_sum = sum(f.^2);
   % f = f ./ sqrt(temp_sum / base_sum);
    
end