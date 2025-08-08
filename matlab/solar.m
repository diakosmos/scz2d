%function crud = solar(stuff)
% boundary conditions: psi = constant (= zero) at theta = 0, pi/2.
%                      d_th omega = zero       at  "  "  "  "  " 
%
%          torque-free:  psi = zero   at   r = rmin, 1
%                  d_r omega = zero    "  "  "  "  "
iterations = 10;

% constants

omega0 = 1;
nu = 1;     % ???/

jmax = 20; % zones in r direction
imax = 30;  % zones in theta direction

dt = 0.001;

rmin = 0.2;  % inner edge of convection zone.

dr = (1-rmin)/(jmax-1);
dth = pi/2 / (imax-1);

kmax = floor(max(imax,jmax/(1-rmin))); % zones in x and y direction for plots -
                               % just needs to give roughly the same or
                               % better resolution than the actual
                               % numerical grid

% r= rmin:(1-rmin)/(imax-1):1;  % is there a better way?
%r1d = rmin + (0:(jmax-1)) * (1-rmin)/(jmax-1);  % safer from round-off?
r1d = linspace(rmin,1,jmax);
th1d = 0:((pi/2)/(imax-1)):pi/2;

r = ones(imax,1) * r1d;
th = th1d' * ones(1,jmax);
[y,x] = pol2cart(th,r);   % this is done "backwards" because it's really
                          % spherical, not cylindrical....
x1d = (0:(kmax-1)) / (kmax - 1);
%xi = ones(kmax,1) * x1d;
%yi = x1d' * ones(1,kmax);
[xi,yi]=meshgrid(x1d,x1d);

% initialize arrays of Omega and Psi

omega = omega0 * ones(imax, jmax);
omega = rand(imax,jmax);
omega = r.^2 .* sin(th);
psi   = zeros(imax, jmax);
%psi=rand(imax,jmax);
s = sin(th);
sinv = [zeros(1,jmax); 1./s(2:imax,:)];
c = cos(th);
lambda = 0.3;
rho = exp(-r / lambda);
d_r_rho = -(1/lambda) * rho;

ell = r .* r .* s .* s .* omega;  % specific angular momentum.
                                  % angular momentum density is ell * rho

%  rho =  ????

while iterations > 0
    iterations = iterations - 1;
% note to self: use "diff" for all of this stuff!!!!    
    for i=(1+1):imax   % "forward centered"
        for j=1:jmax  
            d_th_psi(i,j) = (psi(i,j) - psi(i-1,j))/dth;
            d_th_omega(i,j) = (omega(i,j) - omega(i-1,j))/dth;
            d_th_ell(i,j) = (ell(i,j) - ell(i-1,j))/dth;
        end           
    end
    for j=1:jmax
        d_th_psi(1,j) = d_th_psi(2,j); % lousy
        d_th_omega(1,j) = d_th_omega(2,j);
        d_th_ell(1,j) = d_th_ell(2,2);
    end
    
    for i=1:imax
        for j=(1+1):jmax
            d_r_psi(i,j) = (psi(i,j) - psi(i,j-1))/dr;
            d_r_omega(i,j) = (omega(i,j) - omega(i,j-1))/dr;
            d_r_ell(i,j) = (ell(i,j) - ell(i,j-1))/dr;
        end
    end
    for i=1:imax
        d_r_psi(i,1) = d_r_psi(i,2);
        d_r_omega(i,1) = d_r_omega(i,2);
        d_r_ell(i,1) = d_r_ell(i,2);
    end
    
    for i=1:imax
        for j=2:(jmax-1)
            d_r_r_omega(i,j) = omega(i,j+1) - 2 * omega(i,j) + omega(i,j-1);
            d_r_r_omega(i,j) = d_r_r_omega(i,j) / (dr^2);
        end
    end
    for i=1:imax
        d_r_r_omega(i,1) = d_r_r_omega(i,2);
        d_r_r_omega(i,jmax) = d_r_r_omega(i,jmax-1);
    end

    % find second derivative with respect to theta of omega
    for i=2:(imax-1)
        for j=1:jmax
            d_th_th_omega(i,j) = omega(i+1,j) - 2 * omega(i,j) + omega(i-1,j);
            d_th_th_omega(i,j) = d_th_th_omega(i,j) / (dth^2);
        end
    end
    for j=1:jmax
        d_th_th_omega(1,j) = d_th_th_omega(2,j);
        d_th_th_omega(imax,j) = d_th_th_omega(imax-1,j);
    end
    
    ang_mom_1 = d_th_psi .* d_th_ell;
    ang_mom_2 = d_r_psi  .* d_r_ell;
    ang_mom_3a = nu .* ( d_r_rho .* r .^ 4 .* s .^ 3 .* d_r_omega);
    ang_mom_3b = nu .* rho .* 4 .* r.^3 .* s.^3 .* d_r_omega;
    ang_mom_3c = nu .* rho .* r.^4 .* s.^3 .* d_r_r_omega;
    ang_mom_3 = ang_mom_3a + ang_mom_3b + ang_mom_3c;
    ang_mom_4a = nu .* rho .* r.^2   .*   3 .* s.^2 .* c .* d_th_omega;
    ang_mom_4b = nu .* rho .* r.^2   .*   s.^3 .* d_th_th_omega;
    ang_mom_4 = ang_mom_4a + ang_mom_4b;
    
    %ang_mom_all = ang_mom_1 - ang_mom_2 - ang_mom_3 - ang_mom_4;
    ang_mom_all = -ang_mom_3 - ang_mom_4;
    dt_ang_mom = ang_mom_all ./ (r.^2) .* sinv;
    dtell = dt_ang_mom;
    
    %contour(x,y,ell);
    %axis square;
    %pause;
    %mesh(x,y,ell);
    %axis square;
    %pause;
    elli = griddata(x,y,ell,xi,yi,'cubic');  % not interp2
%    contour(xi,yi,elli);
%    axis square;
%    pause;
    dtelli = griddata(x,y,dtell,xi,yi,'cubic');
    mesh(xi,yi,dtelli);
    %pcolor(dtelli);
    axis square;
    pause;
    
    ell = ell + dtell .* dt;
    omega = ell ./ r ./ r .* sinv .* sinv;   % this really really stinks.
    
    %my_display(x,y,ell);
end
    
% explicit:

%function junk = my_display(x,y,data)
%    datai = griddata(x,y,data,xi,yi,'cubic');
%    mesh(xi,yi,datai);
%    axis square;
%    junk=0;
%end

%end