function null = specsol()

% The program <specsol> is written as a function: null = specsol(),
% but it is best to think of it as a script or program. This
% program calculates the rotation rate of the solar convective zone
% as a function of
% latitude and radius and the coupled meridional (ie poloidal) 
% circulation. It is essentially a 2D hydro code, but the solution
% is by a relaxation method where the nonlinear terms are assumed
% to be in some sense small (XXX check this).
%
% The method of solution is by a pseudospectral
% approach. Discretization is on the spatial grid in the radial
% direction and spectral, on a set of associated Legendre
% functions, in latitude. Azimuthal symmetry is assumed.
%
% This program makes use of the functions chebdif and cheb2bc,
% slightly altered. The original versions are due to
% J.A.C. Weideman and S.C. Reddy, 1998.
%
% (C) 2003 P.T. Williams

%................................................................................
% Boundary Conditions:
%     theta direction:
%         psi           = constant (= zero)  at   theta = (0, pi/2)
%         d_th omega    = zero               at   theta = (0, pi/2)
%         ell, d_th_ell = 0                  at   theta = 0
%         d_th_ell      = 0                  at   theta = pi/2
%                       
%     radial direction:
%         torque-free:
%             psi       = 0                  at   r = (rmin, 1)
%             d_r omega = 0                  at   r = (rmin, 1)
%             d_r ell   = - 2/r ell
%................................................................................
  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  % preliminary initializations

  global omega rmin r r_omega V0 imax jmax nmax c s x y dr s1d c1d ...
      sinv1d leg_n1_si somega ed3

  null = 0;

  % flags:

  visc_f   = 1;                         % flag turn viscosity on
  lambda_f = 1;                         % flag turn lambda effect on
  circ_f   = 0;                         % flag turn meridional circulation on
  varrho_f = 0;                         % flag radially varying density on

  % control parameters:

  jmax = 12;                            % zones in r direction
  imax = 20;                            % zones in theta dithrection
  nmax = 16;                            % number of functions in th-direction
  xmax = 10;                            % number of points in cartesian directions

  dt = 3.0e-6;                          % timestep
  rmin = 0.2;                           % inner edge of convection zone.

  iterations = 1000000;
  plot_step = 100;                      % how often to update plots

  % physical parameters
                                        
  omega0 = 1.0;                         % normalize to this anyway

  if visc_f == 1                        % viscosity, if turned on
    nu = 1.0;                              
  else
    nu = 0.0;
  end

  if lambda_f == 1                      % value of lambda-effect V0
    V0 = 1.0;
  else
    V0 = 0.0;
  end

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  % Stuff to define for the Chebyshev integral:
  % Note that (in TeXspeak) $$ int_{-1}^1 f(x) dx = sum_j( v_j * f(xx_j) ) $$
  % See, eg, Numerical Recipes.

  n = jmax-2;
  xx = mycheb2(n);
  w = pi*(1-xx.^2)/(n+1);  % w, W, v, are for doing integrals - see numrec
  W = sqrt(1-xx.^2);
  v = w./W;

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  % Establish computational mesh.
  %................................................................................
  % First we do the radial mesh. This is done on Chebyshev
  % collocation points:

  % Cheb2bc boundary conditions specification and output format are: 
  % gr = [2 -1 0; 2 -rmin 0];
  % [r1d_us, d2r_us, dr_us, phipr, phimr] = cheb2bc(jmax,gr);
  % Note that vector g = [ap bp cp; am bm cm] as defined on p.25 of Weideman & Reddy
  %
  % Cheb2bc boundary conditions: 
  %     Without Lambda effect:  dr(Omega) = 0
  %     With Lambda effect:     dr(Omega) + V0 * Omega/r = 0
  %
  % But beware, this is a bit tricky, because cheb2bc assumes boundaries are
  % at +-1, so you have to fool it! In particular, this will NOT work:
  % gr = [(V0) 1 0; (V0) rmin 0];
  %
  % The trick is that not only do you have to rescale the spatial grid & the derivative
  % matrices, but you also have to rescale the boundary conditions.
  %
  % Note that "_us" suffix on variable names indicate an *unscaled*
  % variable, ie one that has not been scaled to the [rmin,1]
  % interval.
                                      
  rmult = 0.5 * (1.0 - rmin);           % scale from [1,-1] to [rmin,1]
  gr = [(V0)*rmult 1 0; (V0)*rmult rmin 0];         
  [r1d_us, d2r_us, dr_us, phipr, phimr] = cheb2bc(jmax,gr);

  %[r1d_us, dnr_us] = chebdif(jmax,2);  % bc appropriate for ell?
  %dr_us = dnr_us(:,:,1);
  %d2r_us = dnr_us(:,:,1);

  r1d = rmin + rmult * (r1d_us + 1.0);  % r1d is 1D (vector) radius mesh
  dr = dr_us ./ rmult;                  % dr is radial derivative matrix
  d2r = d2r_us ./ rmult^2;              % d2r is 2nd radial deriv matrix

  sr = ones(nmax,1) * r1d';             % sr is Spectral Radius matrix
  srinv = 1./sr;                        % srinv is XXX NASTY. Hope not to use.    

  % BTW, uncomment the next three lines if you want to check for yourself that you
  % are taking the right derivatives etc. (Take a function we know
  % the derivative of, and compare with derivative given by dr matrix.):
  %s = sin(th);
  %dths = dth * s;
  %drs = (dr * s')';

  %................................................................................
  % Next step is building grid in theta-direction, but this is potentially more
  % subtle since in fact we do an expansion in associated Legendre functions.
  %
  % The only real use for such a spatial grid is to project initial
  % conditions from spatial to spectral, if they are specified in
  % the spatial domain, and to back-transform final results onto
  % spatial mesh to create nice plots. (This is trivial of course.)
  %
  % Spatial theta grid should really be done using Gauss-Legendre
  % etc, but this is not too important if the algorithm just
  % transforms once each way at the beginning and end of the
  % entire run, as described...... 
  %
  % So here we stick with a uniform mesh.

  th1d = linspace(0,pi,imax);

  %................................................................................
  % We are now ready to create the full spatial mesh and associated functions.

  [r,th] = meshgrid(r1d,th1d);

  % Cartesian meshes, to interpolate stuff onto a Cartesian grid if desired:
                                        
  [y,x] = pol2cart(th,r);               % this is done "backwards" because it's really
                                        % spherical, not cylindrical....
  x1d = linspace(0,1,xmax);
  y1d = linspace(-1,1,2*xmax);
  [xi,yi] = meshgrid(x1d,y1d);

  % various associated functions (eg sin(theta), etc)

  s   = sin(th.*1.0);
  s1d = sin(th1d .* 1.0);
  c   = cos(th);
  c1d = cos(th1d .* 1.0);

  % 1/r and 1/theta meshes, including XXX NASTY sinv
  % Note that sinv is formally infinite at th = 0, pi, but set to zero here.

  rinv = 1./r;
  sinv = [zeros(1,jmax); 1./s(2:imax-1,:); zeros(1,jmax)]; 
  sinv1d = [0, 1./s1d(2:imax-1), 0];

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  % Now, initialize physical variable (eg omega(i,j)) meshes:

  initialize(0);                        % initialize omega

  % initialize density 1D matrix

  rho1d = ones(jmax,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop:
  mod = 0;   % "mod" for "modulo": ie, a way to plot only, say, every 100th iteration.
  step = 0;
  while step < iterations
    step = step+ 1;

    % back transform:
    for j=1:jmax,
      omega(:,j) = leg_n1_si' * somega(:,j);
    end

    % plot:
    if mod==0
      mod = plot_step;
      mesh(x,y,omega);
      pause(0.01);
    end
    mod = mod-1;

    % Plotting may be done in cartesian also:    
    %    omegacart = griddata(x,y,omega,xi,yi,'cubic');
    %    mesh(xi,yi,omegacart);
    %    pause(0.01);

    dr_somega = (dr * somega')';			 
    %term_a = srinv.^4 .* (dr *  (sr.^4 .* dr_somega)')';
    visc_a = nu * ((d2r * somega')'  + 4 * srinv .* dr_somega); 
    % term_a = 0;
    %			 mesh(somega); pause;
    %			 mesh(term_a); pause;
    visc_b = nu * (-srinv.^2 .* (ed3 .* somega));
   
    lamb_a = nu * (dr * ( (somega .* srinv) * V0)')'   +   4 * nu * (srinv.^2 .* somega) * V0;


    somega = somega + dt *  ( visc_a + visc_b  + lamb_a);
    
    % checking conservation of angular momentum:

    %L1 = rho .* x.^3 .* 
    om1 = somega(1,:)';   % omega_1 - ie in spectral space - as a function of r.
    om3 = somega(3,:)';
    f1 =  (rho1d .* r1d.^3 .* om1);
    f2 =  (rho1d .* r1d.^3 .* om3);
    int1 = 2*pi * (-16/15) * sum( v .* f1(2:jmax-1));
    int2 = 2*pi * (16/35) * sum(v .* f2(2:jmax-1));
    if step == 1 
      Ltot = int1 + int2;
    end
    Ltotnew = int1 + int2;
    Lhist(step) = Ltotnew;
    %ratio = Ltot / Ltotnew
    %somega = somega * (Ltot / Ltotnew);
   
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function dum1 = initialize(dum2)
  dum1 = 0
  global omega rmin r r_omega V0 imax jmax nmax c s x y dr s1d c1d ...
      sinv1d leg_n1_si somega ed3
  
  % The solution should - we hope - asymptote to a unique solution,
  % modulo the conserved quantities, eg total angular momentum and
  % energy (XXX check - how does energy play into this?). If this
  % is true, then initial conditions are largely irrelevant.
  %
  % The main consideration is that the initial condition is
  % consistent with the boundary conditions. Otherwise, due to the
  % use of the cheb2bc routine, we may see gross errors in initial
  % behavior, which may seed numerical instabilities.

  %................................................................................
  % Radial dependence:
  % Start with 3rd-order poly that has zero derivative at rmin and 1.0.

  a=rmin; b=1.0;
  r_omega = 100 * (2*r.^3 - 3*(a+b)*r.^2 + 6 * a * b * r - 1.4);    

  % Now fix radial dependence of omega so that boundary conditions are satisfied
  % even if lambda term is included:

  fi =  (0.5) * (r - a).^2 / (b - a);  % has derivative 1 at r_max, zero derivative at rmin, also equal to zero ar rmin.
  fo = -(0.5) * (r - b).^2 / (b - a);  % has derivative 1 at rmin, zero at rmax, also equal to zero at rmax.
                                       %r_omega = r_omega - (V0 * r_omega(1,1) / b) * fo - (V0 * r_omega(1,jmax) / a) * fi;
  alpha = - (r_omega(1,1)/b) / (1 + (0.5*(b-a)/b));
  beta  = - (r_omega(1,jmax)/a) / (1 - (0.5 * (b-a)/a));
  r_omega = r_omega + alpha * fi + beta * fo;

  omega = r_omega .* (c + s);
  omega = omega * (-1);
  %ell = omega .* r.^2 .* s.^2;
  %ellibc = ell(:,imax);

  mesh(x,y,(dr * omega')' + V0 * omega./r);

  % mesh(x,y,omega);
  pause;

  % Project initial conditions onto associated Legendre funcions

  omega_start = omega;

  % Note that if g(th) = sum( g_n * P_n^1(cos(th)) / sin(th)
  % then
  % g_n = (2n+1)/(2n * (n+1)) * int_0^pi g(th) P_n^1(cos(th)) * sin^2(th) d(th)

  % Integral best done using Gauss-Legendre, but here is a cheap way:

  nn = 1:nmax;
  norm = (2 * nn + 1) ./ ((2*nn).*(nn+1)); % normalization for P_n^1.
  ed3_vec = (nn-1) .* (nn+2);  % e-value of operator (1/s^3) * d_th( s^3 * d_th( P_n^1/s )))
  ed3 = ed3_vec' * ones(1,jmax);

  for n= 1:nmax,
    leg_nm = legendre(n,c1d);
    leg_n1(n,:) = leg_nm(2,:);
  end
  for n=1:nmax, 
    for i=2:imax-1,
      leg_n1_si(n,i) = leg_n1(n,i) * sinv1d(i); 
    end
    leg_n1_si(n,1) = -n * (n+1) / 2;
    leg_n1_si(n,imax) = (-1)^n * n * (n+1)/2;
  end
  somega = zeros(nmax,jmax);
  dth = pi/(imax-1);
  for n=1:nmax,
    %   for j=1:jmax,  % loop over r
    %   Note that P_1^1(cos th) / sin( th) = -1, so this takes care of constant.
    somega(n,:) = norm(n) * (leg_n1(n,:) .* s1d.^2) * omega  * dth;
    %      somega(n,j) = norm(n) * (leg_n1(n,:) .* s1d.^2) * omega(:,j);
    %   end			   
  end
