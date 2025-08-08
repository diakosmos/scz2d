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
clear;
iterations = 1000;

% constants

omega0 = 1;
nu = 1;     % ???/

n = 20; % zones in r direction
m = 30;  % zones in theta dithrection

dt = 0.01;

rmin = 0.2;  % inner edge of convection zone.

%kmax = floor(max(imax,jmax/(1-rmin))); % zones in x and y direction for plots -
                               % just needs to give roughly the same or
                               % better resolution than the actual
                               % numerical grid

%[r1d_us,dr_us] = chebdif(jmax,2);
gr=[0,1,0;0,1,0];
[r1d_us,d2r_us,dr_us,omega_outer,omega_inner]=cheb2bc(n,gr);
% r1d_us ranges from -1 to 1. rescale to get r. remember this when you take
% derivatives.
%[th1d_us,dth_us] = chebdif(imax,2);
gth=[0,1,0;0,1,0];
[th1d_us,d2th_us,dth_us,omega_thpi,omega_th0]=cheb2bc(m,gth);

r1d = r1d_us .* ((1.0 - rmin)/2)  +  (1.0 + rmin)/2;
%dr = dr_us(:,:,1) ./ ((1.0 - rmin)/2);
dr = dr_us ./ ((1.0-rmin)/2);
%d2r = dr_us(:,:,2) ./ ((1.0 - rmin)/2)^2;
d2r = d2r_us ./ ((1.0-rmin)/2)^2;

th1d = (th1d_us + 1)./2  .* pi;
%dth = dth_us(:,:,1) ./ (pi / 2);
dth = dth_us ./ (pi/2);
%d2th = dth_us(:,:,2) ./ (pi / 2)^2;
d2th = d2th_us ./ (pi/2)^2;

[r,th] = meshgrid(r1d,th1d);
%[th,r] = meshgrid(th1d,r1d);
dth_flat=mymult(eye(n),dth);
dr_flat=mymult(dr,eye(m));
d2th_flat=mymult(eye(n),d2th);
d2r_flat=mymult(d2r,eye(m));
th_flat = flatten(th);
r_flat = flatten(r);
[y,x] = pol2cart(th,r);   % this is done "backwards" because it's really
                          % spherical, not cylindrical....
%x1d = linspace(0,1,kmax);
%y1d = linspace(-1,1,2*kmax);
%[xi,yi]=meshgrid(x1d,y1d);

s = sin(th.*1.0);
c= cos(th);
sinv = [zeros(1,n); 1./s(2:m-1,:); zeros(1,n)];

p5 = legendre(5,c);
p51 = squeeze(p5(2,:,:));
p51_os = p51 .* sinv;


% initialize arrays of Omega and Psi

% Uncomment the next three lines if you want to check for yourself that you
% are taking the right derivatives etc.:
%s = sin(th);
%dths = dth * s;
%drs = (dr * s')';

a=rmin; b=1.0;
%omega = r .^ 2 .* cos(th); % initialize with some nontrivial, uh, r- and th- dependence.
r_omega = 1.0 * (2*r.^3 - 3*(a+b)*r.^2 + 6 * a * b * r + 0.0);
d2r_r_omega = 12* r - 6*(a+b);
% 3rd-order poly that has zero derivative at rmin and 1.0.
omega = r_omega .* c;
d2r_omega = d2r_r_omega .* c;

d2r=popup(d2r_flat * flatten(omega),m,n);



% find angular momentum
    
    ell = r.^2 .* s.^2 .* omega;

% for now keep rho = nu = 1
r_f=flatten(r);
s_f=flatten(s);
c_f = flatten(c);
sinv_f=flatten(sinv);
    ang_mom_3a = diag(4 .* r_f .* s_f.^2) * dr_flat;
    ang_mom_3b = diag(r_f.^2 .* s_f.^2) * d2r_flat;
ang_mom_4a = 0;
      ang_mom_4a = diag(3 .* s_f .* c_f) * dth_flat;
    ang_mom_4b = diag(s_f.^2) * d2th_flat; 

    dt_ell = ang_mom_3a + ang_mom_3b + ang_mom_4a + ang_mom_4b;

    dt_omega = diag(r_f .\ r_f .\  sinv_f ) * dt_ell;  % that nasty sinv is the root of all evil!
    %dt_omega = d2r_flat + d2th_flat;
    f=omega;
    ff=flatten(f);
    g=ell;
    gg=flatten(g);
dt = 0.01;

mymat=eye(n*m)-dt*dt_omega;
mymat2=eye(n*m)-dt*dt_ell;
while iterations > 0
    iterations = iterations - 1;

%    mesh(f);
mesh(x,y,g);
pause(0.1);
    gg = mymat2 \ gg;
    g=popup(gg,m,n);
    
    
    % this next step is the way you would do it explicitly..
    %omega = omega + dt_omega .* dt;
    
    %mymat = eye(
    
    %elli = griddata(x,y,ell,xi,yi,'cubic');
    %mesh(xi,yi,elli);

% update angular momentum    
     
%    ell = ell + dt_ell * dt;

% invert to get omega -- note that omega at theta=0 is undefined here, but irrelevant as well.    
    
end