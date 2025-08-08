%function crud = solar(stuff)
% boundary conditions: psi = constant (= zero) at theta = 0, pi/2.
%                      d_th omega = zero       at  "  "  "  "  " 
%                      ell, d_th_ell = 0 at    theta = 0
%                           d_th_ell = 0 at    theta = pi/2
%                       
%
%          torque-free:  psi = zero   at   r = rmin, 1
%                  d_r omega = zero    "  "  "  "  "
%                  d_r ell = - 2/r ell

iterations = 10;

% constants

omega0 = 1;
nu = 1;     % ???/

jmax = 15; % zones in r direction
imax = 15;  % zones in theta dithrection

dt = 0.001;

rmin = 0.2;  % inner edge of convection zone.

dr = (1-rmin)/(jmax-1);
dth = pi/2 / (imax-1);

kmax = floor(max(imax,jmax/(1-rmin))); % zones in x and y direction for plots -
                               % just needs to give roughly the same or
                               % better resolution than the actual
                               % numerical grid

r1d = linspace(rmin,1,jmax);
th1d = linspace(0,pi/2,imax);
[r,th]=meshgrid(r1d,th1d);
[y,x] = pol2cart(th,r);   % this is done "backwards" because it's really
                          % spherical, not cylindrical....
x1d = linspace(0,1,kmax);
[xi,yi]=meshgrid(x1d,x1d);

% initialize arrays of Omega and Psi

omega = omega0 * ones(imax, jmax);
omega = rand(imax,jmax);
omega = r.^2 .* sin(th);
psi   = zeros(imax, jmax);
%psi=rand(imax,jmax);
s = sin(th);
sinv = [zeros(1,jmax); 1./s(2:imax,:)];  % sinv is defined here to be zero when sin(th) = 0
c = cos(th);
rinv = 1 ./ r;
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
        d_th_ell(1,j) = 0; %d_th_ell(2,2);
        d_th_ell(imax,j) = 0;    
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
        d_r_ell(i,1) = -2 .* rinv(i,1) .* ell(i,1);
        d_r_ell(i,jmax) = -2 * rinv(i,jmax) * ell(i,jmax);
    end
    
    for i=1:imax
        for j=2:(jmax-1)
            d_r_r_omega(i,j) = omega(i,j+1) - 2 * omega(i,j) + omega(i,j-1);
            d_r_r_omega(i,j) = d_r_r_omega(i,j) / (dr^2);
            
            d_r_r_ell(i,j) = ell(i,j+1) - 2 * ell(i,j) + ell(i,j-1);
            d_r_r_ell(i,j) = d_r_r_ell(i,j) / (dr^2);
        end
    end
    for i=1:imax
        d_r_r_omega(i,1) = d_r_r_omega(i,2);
        d_r_r_omega(i,jmax) = d_r_r_omega(i,jmax-1);
        
        d_r_r_ell(i,1) = d_r_r_ell(i,2);
        d_r_r_ell(i,jmax) = d_r_r_ell(i,jmax-1);
    end

    % find second derivative with respect to theta of omega
    for i=2:(imax-1)
        for j=1:jmax
            d_th_th_omega(i,j) = omega(i+1,j) - 2 * omega(i,j) + omega(i-1,j);
            d_th_th_omega(i,j) = d_th_th_omega(i,j) / (dth^2);
            d_th_th_ell(i,j) = (ell(i+1,j) - 2 * ell(i,j) + ell(i-1,j)) ./ (dth^2);
        end
    end
    for j=1:jmax
        d_th_th_omega(1,j) = d_th_th_omega(2,j);
        d_th_th_omega(imax,j) = d_th_th_omega(imax-1,j);
        d_th_th_ell(1,j) = d_th_th_ell(2,j);
        d_th_th_ell(imax,j) = d_th_th_ell(imax-1,j);
    end

    ang_mom_3a = d_r_r_ell - 2 .* rinv .* rinv .* ell;
    ang_mom_3b = (d_r_ell - 2 .* rinv .* ell)./lambda;
    ang_mom_4a = -rinv .* rinv .* ( d_th_th_ell - c .* sinv .* d_th_ell - s .* ell);
    
    ang_mom_all = - rho .* nu .* (ang_mom_3a + ang_mom_3b + ang_mom_4a);

    dtell = ang_mom_all;
    
    % enforce BC:
        for i=2:(imax-1)
        for j=1:jmax
            d_th_th_omega(i,j) = omega(i+1,j) - 2 * omega(i,j) + omega(i-1,j);
            d_th_th_omega(i,j) = d_th_th_omega(i,j) / (dth^2);
            d_th_th_ell(i,j) = (ell(i+1,j) - 2 * ell(i,j) + ell(i-1,j)) ./ (dth^2);
        end
    end
    for j=1:jmax
        d_th_th_omega(1,j) = d_th_th_omega(2,j);
        d_th_th_omega(imax,j) = d_th_th_omega(imax-1,j);
        d_th_th_ell(1,j) = d_th_th_ell(2,j);
        d_th_th_ell(imax,j) = d_th_th_ell(imax-1,j);
    end

%end
     pause;
    %pause;
%    axis square;
%   dtelli = griddata(x,y
    elli = griddata(x,y,ell,xi,yi,'cubic');  % not interp2
%     mesh(xi,yi,elli);
   %contour(x,y,ell);
     %pcolor(dtelli
   contour(xi,yi,elli);,dtell,xi,yi,'cubic');
    %axis square;
    %axis );
    axis square;
    pause;
    
    ell = ell - dtell .* dt;
    ell(1,:) = 0;
%    omega = ell ./ r ./ r .* sinv .* sinv;   % this really really stinks.
       %pause;
 
   %mesh(x,y,ell);
     %my_display(x,y,ell);
end
    
% explicit:

%function junk = my_display(x,y,data)
%    datai = griddata(x,y,data,xi,yi,'cubic');
%    mesh(xi,yi,datai);
%    axis square;
%    junk=0;
%end
