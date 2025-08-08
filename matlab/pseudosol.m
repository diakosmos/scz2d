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

jmax = 40; % zones in r direction
imax = 100;  % zones in theta dithrection

dt = 0.001;

rmin = 0.2;  % inner edge of convection zone.

%kmax = floor(max(imax,jmax/(1-rmin))); % zones in x and y direction for plots -
                               % just needs to give roughly the same or
                               % better resolution than the actual
                               % numerical grid

[r1d_us,dr_us] = chebdif(jmax,2);
% r1d_us ranges from -1 to 1. rescale to get r. remember this when you take
% derivatives.
[th1d_us,dth_us] = chebdif(imax,2);

r1d = r1d_us .* ((1.0 - rmin)/2)  +  (1.0 + rmin)/2;
dr = dr_us(:,:,1) ./ ((1.0 - rmin)/2);
d2r = dr_us(:,:,2) ./ ((1.0 - rmin)/2)^2;
th1d = (th1d_us + 1)./2  .* pi;
dth = dth_us(:,:,1) ./ (pi / 2);
d2th = dth_us(:,:,2) ./ (pi / 2)^2;

[r,th] = meshgrid(r1d,th1d);
%[th,r] = meshgrid(th1d,r1d);

[y,x] = pol2cart(th,r);   % this is done "backwards" because it's really
                          % spherical, not cylindrical....
%x1d = linspace(0,1,kmax);
%y1d = linspace(-1,1,2*kmax);
%[xi,yi]=meshgrid(x1d,y1d);

s = sin(th.*1.0);
c= cos(th);
sinv = [zeros(1,jmax); 1./s(2:imax-1,:); zeros(1,jmax)];

% initialize arrays of Omega and Psi

% Uncomment the next three lines if you want to check for yourself that you
% are taking the right derivatives etc.:
%s = sin(th);
%dths = dth * s;
%drs = (dr * s')';

omega = r .^ 2 .* c; % initialize with some nontrivial, uh, r- and th- dependence.
a=rmin; b=1.0;
r_omega = 1.0 * (2*r.^3 - 3*(a+b)*r.^2 + 6 * a * b * r + 0.0);    
% 3rd-order poly that has zero derivative at rmin and 1.0.
omega = r_omega .* c;

% main loop:

while iterations > 0
    iterations = iterations - 1;

    mesh(x,y,omega);

    pause;
    
    dr_omega = (dr * omega')';
    dth_omega = dth * omega;
    d2r_omega = (d2r * omega')';
    d2th_omega = d2th * omega;
    d2rth_omega = (dr * (dth * omega)')';
    d2thr_omega = dth * (dr * omega')';

% enforce Dirichlet boundary conditions:

    dr_omega(:,1) = 0;
    dr_omega(:,jmax) = 0;
    dth_omega(1,:) = 0;
    dth_omega(imax,:)=0;

% find angular momentum
    
    ell = r.^2 .* s.^2 .* omega;

% for now keep rho = nu = 1

    ang_mom_3a = 4 .* r .* s.^2 .* dr_omega;
    ang_mom_3b = r.^2 .* s.^2 .* d2r_omega;
    ang_mom_4a = 3 .* s.^2 .* sinv .* c .* dth_omega;
    ang_mom_4b = s.^2 .* d2th_omega; 

    dt_ell = ang_mom_3a + ang_mom_3b + ang_mom_4a + ang_mom_4b;

    dt_omega = dt_ell ./ r ./ r .* sinv .* sinv;
    
    omega = omega + dt_omega .* dt;
    
    %elli = griddata(x,y,ell,xi,yi,'cubic');
    %mesh(xi,yi,elli);

% update angular momentum    
     
%    ell = ell + dt_ell * dt;

% invert to get omega -- note that omega at theta=0 is undefined here, but irrelevant as well.    
    
end