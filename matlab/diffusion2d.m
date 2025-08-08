%clear;
n=30; m=50;
gx=[0,1,0;0,1,0];
[x,d2x,d1x,xm,xp]=cheb2bc(n,gx);
gy=[0,1,0;0,1,0];
[y,d2y,d1y,ym,yp]=cheb2bc(m,gy);
%[x,d] = chebdif(n,2);
%d1x = d(:,:,1); d2x = d(:,:,2);
%[y,d] = chebdif(m,2);
%d1y = d(:,:,1); d2y = d(:,:,2);
dxxx=mymult(d1x,eye(m));
dyyy=mymult(eye(n),d1y);
d2xxx=mymult(d2x,eye(m));
d2yyy=mymult(eye(n),d2y);

[xx,yy]=meshgrid(x,y);
xxx=flatten(xx);
yyy=flatten(yy);
ax = 0.4;  % +- zero of derivative of function below - other zeros of derivative are +-1.
fx= 3*xx.^5 - 5*(ax^2 + 1)*xx.^3 + 15 * ax^2 .* xx;  % polynomial that is symmetric, nontrivial, and
       % satisfies boundary conditions
fx = fx + (1 - xx.^2);   % make its integral nonzero, to test conservation, since integral of f should
% be conserved, with bc that f's derivative is zero on boundaries.
%fx = xx;
ay = 0.2;
fy = 3*yy.^5 - 5*(ay^2 +1)*yy.^3 + 15 * (ay^2) * yy;
fy = fy + (1 - yy.^4);
%fy = yy.^2;
f = fx + fy;
dxf=popup(dxxx * flatten(f),m,n);
dyf=popup(dyyy * flatten(f),m,n);
d2xf=popup(d2xxx * flatten(f),m,n);
d2yf=popup(d2yyy * flatten(f),m,n);
pause;

%base_sum_1 = sum(sum(f));       
%base_sum_2 = sum(sum(f.^2));
%hold off;
%plot(x,f);       
%    hold on;   
ff=flatten(f);
d2 = d2xxx + d2yyy;
dt = 0.01;

mymat=eye(n*m)-dt*d2;
%mymatinv = mymat \ 1.0;
iterations = 1000; i=0;
while i < iterations
    i = i + 1;
 %   base_sum(i) = sum(f);
 %   f = f ./ (base_sum(i) / base_sum_1);
 %   plot(x,f);
    
    %pause;
    
    %d2f = d2 * f;

    %dtf = d2f * dt;

    %f = f + dtf;
    
    %mymat = eye(n) - dt * d2;
    %invmymat = 1/mymat;
    %f = mymat \ f;
    mesh(xx,yy,f);
    pause(0.1);
    ff = mymat \ ff;
    f=popup(ff,m,n);
   % temp_sum = sum(f.^2);
   % f = f ./ sqrt(temp_sum / base_sum);
    
end