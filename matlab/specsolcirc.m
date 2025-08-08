function null = specsolcirc()
%
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
% direction and spectral, on a set of (modified) associated Legendre
% functions, in latitude. Azimuthal symmetry is assumed.
%
% There are two functions that specify the flow: Omega and Psi, ie
% the angular velocity and the streamfunction. Latitudinally, Omega
% is expressed spectrally in terms of functions f_n, and Psi is
% expanded in functions g_n.
%
% Our particular choice for f_n and g_n is P_n^1(c)/s and
% P_n^1(c)*s, respectively, where c=cos(theta) and s=sin(theta).
%
% The radial grid is chosen to lie on Chebyshev collocation points.
%
% For radial derivatives on the spatial grid,
% this program makes use of the functions chebdif and cheb2bc,
% slightly altered. The original versions are due to
% J.A.C. Weideman and S.C. Reddy, 1998.
%
% (C) 2003/2004 P.T. Williams
%
% For further information, contact the author at ptw@lanl.gov or
% petwil@astro.as.utexas.edu, or search astro-ph for publications
% by the author related to this work.

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
%             d_r omega + (V0/r)*Omega = 0   at   r = (rmin, 1)
%             d_r ell   = - 2/r ell  XXX is this right?
%................................................................................

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  % preliminary initializations

  null = 0.0;                           % function must return a value
  global imax jmax nmax norm leg_n1 leg_n1_si s1d

  % flags:

  visc_f   = 1;                         % flag to turn viscosity on
  lambda_f = 1;                         % flag to turn lambda effect on
  circ_f   = 0;                         % flag to turn meridional circulation on
  varrho_f = 0;                         % flag to turn radially varying density on

  % control parameters:

  jmax = 12;                            % zones in r direction
  imax = 20;                            % zones in theta dithrection
  nmax = 16;                            % number of functions in th-direction
  xmax = 10;                            % number of points in cartesian directions

  dt = 3.0e-6;                          % timestep
  rmin = 0.2;                           % inner edge of convection zone.

  iterations = 1000000;
  plot_step = 1000;                      % how often to update plots

  % physical parameters
                                        
  omega0 = 1.0;                         % normalize to this anyway

  if visc_f == 1                        % viscosity, if turned on
    nu = 1.0;                              
  else
    nu = 0.0;
  end

  if lambda_f == 1                      % lambda-effect V0,
    V0 = -1.0;                           % if turned on
  else
    V0 = 0.0;
  end

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  % Stuff to define for the Chebyshev integral (for the radial direction):
  % Note that (in TeXspeak) $ int_{-1}^1 f(x) dx = sum_j( v_j * f(xx_j) ) $
  % See, eg, Numerical Recipes.

  n = jmax-2;
  xx = mycheb2(n);
  w = pi*(1-xx.^2)/(n+1);  % w, W, v, are for doing integrals - see numrec
  W = sqrt(1-xx.^2);
  v = w./W;
  clear n w W

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  % Establish computational mesh. There are two: omega and psi.
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
  srinv = 1./sr;                        % srinv is just 1/r; ok
                                        % since rmin > 0.    

  % BTW, mess with the next three lines if you want to check for yourself that you
  % are taking the right derivatives etc. (Take a function we know
  % the derivative of, and compare with derivative given by dr matrix.):
  %s = sin(th);
  %dths = dth * s;
  %drs = (dr * s')';

  %................................................................................
  % Next step is building grid in theta-direction. This is potentially
  % subtle since in fact we do an expansion in associated Legendre functions.
  % (sort of)
  %
  % The only real use for such a spatial grid is to project initial
  % conditions from spatial to spectral, if they are specified in
  % the spatial domain, and to back-transform final results onto
  % spatial mesh to create nice plots. (This final step is trivial of course.)
  %
  % Spatial theta grid should really be done using Gauss-Legendre
  % etc, but this is not too important if the algorithm just
  % transforms once each way at the beginning and end of the
  % entire run, or even better, starts out in spectral and
  % transforms only once, to spatial, at the end, as described...... 
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
  c   = cos(th.*1.0);
  c1d = cos(th1d .* 1.0);

  % 1/r mesh, and XXX NASTY sinv (1 / sin(theta)) mesh.
  % Note that sinv is formally infinite at th = 0, pi, but set to zero here.
  % If sinv is actually used anywhere, it is a red flag that the
  % numerical method needs some work.

  rinv = 1./r;
  sinv = [zeros(1,jmax); 1./s(2:imax-1,:); zeros(1,jmax)]; 
  sinv1d = [0, 1./s1d(2:imax-1), 0];

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  % Legendre functions.
  % Note that legendre(n,x) returns a column vector of the
  % associated legendre functions of degree n and order m = 0..n,
  % evaluated at x (x = cos(theta), typically).

  %................................................................................  
  % Create leg_n1(1:nmax,1:imax): This contains legendre
  % P_n^1(cos(th)) for all n up to nmax, for all points on spatial
  % theta grid.

  for n= 1:nmax,
    leg_nm = legendre(n,c1d);           % Recall c1d is cos(th) on 1d theta-mesh.
    leg_n1(n,:) = leg_nm(2,:);          % Pick out 1st order - ie P_n^1(x).
  end

  %................................................................................
  % Create P_n^1/s. This is a little tricky at (theta=0,pi), since we
  % get a 0/0 problem there, so we have to insert by hand the
  % asymptotic value.

  for n=1:nmax, 
    for i=2:imax-1,
      leg_n1_si(n,i) = leg_n1(n,i) * sinv1d(i); 
    end
    leg_n1_si(n,1) = -n * (n+1) / 2;
    leg_n1_si(n,imax) = (-1)^n * n * (n+1)/2;
  end

  %................................................................................
  % some other ancillary stuff related to normalization

  nn = 1:nmax;
  norm = (2 * nn + 1) ./ ((2*nn).*(nn+1)); % normalization for P_n^1.
  ed3_vec = (nn-1) .* (nn+2);  % e-value of operator (1/s^3) * d_th( s^3 * d_th( P_n^1/s )))
  ed3 = ed3_vec' * ones(1,jmax);

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  % Now, initialize physical variable (eg omega(i,j)) meshes:

%  initialize(0);                        % initialize omega
%  psi = omega * 0.0;                    % initialize psi - initially zero.

  % initialize density 1D matrix

  rho1d = ones(jmax,1);
  %................................................................................
  % Radial dependence:
  % Start with r.e. method:
  % 3rd-order poly that has zero derivative at rmin and rmax.
  % Actually you can put anything here so long as it has zero
  % derivative at rmin and rmax.

  a=rmin; b=1.0;
  r_omega_0 = -1.0 * (2*r1d.^3 - 3*(a+b)*r1d.^2 + 6 * a * b * r1d - 1.4);    
  %r_omega_0 = r1d * 0.0 + 1;

  % Now fix radial dependence of omega so that boundary conditions are satisfied
  % even if lambda term is included.
  %
  % Use functions fi and fo defined below, on interval [a,b] (where
  % a = rmin, b = rmax = 1.0):

  fi =  (0.5) * (r1d - a).^2 / (b - a);  % fi is 0 at a, 0.5(b-a) at b.
                                       % derivative of fi is 0
                                       % at a, 1 at b.   
  fo = -(0.5) * (r1d - b).^2 / (b - a);  % fo is -0.5(b-a) at a, 0 at b.
                                       % derivative of fo is 1 at a,
                                       % 0 at b.
  % Now given the initial function r_omega_0, which has zero
  % derivative on the boundaries [a,b], and the functions fi and fo,
  % there exist the numbers alpha and beta so that the radial
  % function $r_omega = r_omega_0 + alpha * fi + beta * fo$ will
  % satisfy the boundary condition that \p_r Omega + (V0/r) Omega =
  % 0 at rmin and rmax (ie a and b). This works, as you can verify
  % with some only slightly tedious algebra:                    

  alpha = - (r_omega_0(1)/b) / (1/V0 + (0.5*(b-a)/b));
  beta  = - (r_omega_0(jmax)/a) / (1/V0 - (0.5 * (b-a)/a));
  r_omega = r_omega_0 + alpha * fi + beta * fo

  %................................................................................
  % Theta-dependence:
  % Previously I used some randomish function like cos(th) +
  % sin(th), and then projected this onto the spectral
  % space. However, it is easier to just start off with the
  % spectral form, and so that is what I intend to do now. 
                                        
  s_th_omega = sqrt(norm)' .* (1.0 * 2 * (rand([nmax,1]) - 0.5));   % random spectral
                                            % initialization
  s_th_omega(1) = -1.0;                     % except for 1 and 3,
  s_th_omega(3) = 1.0;                      % so that total ang mom
                                            % is positive.
  %................................................................................
  % Combine by just multiplying. No need to get any more fancy than
  % this.
  
  somega = s_th_omega * r_omega';

  % Note that if g(th) = sum( g_n * P_n^1(cos(th)) / sin(th)
  % then
  % g_n = (2n+1)/(2n * (n+1)) * int_0^pi g(th) P_n^1(cos(th)) * sin^2(th) d(th)

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  omega = om_b_trans(somega);

  %................................................................................
  % Check on satisfying boundary conditions &c.
  % This should produce something smooth looking. If it looks
  % jaggedy, it is due to artifacts in cheb2bc, and this indicates
  % that the bounardy conditions are messed up somewhere:
  mesh(x,y,(dr * omega')' + V0 * omega./r);
  pause;

  scaler = 1./norm' * ones(1,jmax);  % created solely for
                                     % displaying spectrum - it's
                                     % better to plot strengths of
                                     % normalized spectrum than unnormalized.
  evenodd = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop:
  mod = 0;   % "mod" for "modulo": ie, a way to plot only, say, every 100th iteration.
  step = 0;
  while step < iterations
    step = step+ 1;

    % back transform:
    omega = om_b_trans(somega);
%    for j=1:jmax,
%      omega(:,j) = leg_n1_si' * somega(:,j);
%    end

    % plot:
    if mod==0
      mod = plot_step;
      evenodd = evenodd * -1;
      if (evenodd > 0)
 %       pcolor(scaler .* somega);
 %       colorbar('vert');
      else
%        mesh(x,y,omega);
%      contour(x,y,omega);
    % Plotting may be done in cartesian also:    
%        omegacart = griddata(x,y,omega,xi,yi,'cubic');
%        mesh(xi,yi,omegacart);
      end
      pause(0.01);
    end
    mod = mod-1;

    dr_somega = (dr * somega')';			 

    %................................................................................
    % viscous terms - XXX note that if rho is not constant, there
    % are additional terms that I have not included yet !!!!
    %
    % (1\over \tilde \nu) \partial_\tau\Omega_n(x) =
    %        {1\over x^4} \partial_x(x^4 \partial_x \Omega_n(x)) -
    %      - {1\over x^2} (n-1)(n+2) \Omega_n(x)
    %
    % Radial derivative term:

    visc_a = nu * ((d2r * somega')'  + 4 * srinv .* dr_somega); 

    % Latitudinal derivative term:

    visc_b = nu * (-srinv.^2 .* (ed3 .* somega));
   
    %................................................................................
    % lambda-effect terms
    %
    % {1\over x^4 s^2} \partial_x( x^4 s^2 \rho \tilde\nu {\Omega
    % \over x} V0 =
    % \partial_x(\rho \tilde\nu {\Omega \over x} V0) + (4/x) \rho
    % \tilde\nu (\Omega / x) V0

    lamb_a = nu * (dr * ( (somega .* srinv) * V0)')'   +   4 * nu * (srinv.^2 .* somega) * V0;

    %................................................................................
    % circulation (ie streamfunction, psi) terms
    %
    %  \partial_t(\rho \Omega) = 
    %       { 1 \over r^4 s^4} (\partial_\theta \psi) \partial_r(r^2 s^2 \Omega)   -
    %       { 1 \over r^4 s^3} (\partial_r \psi)  \partial_\theta(r^2 s^2 \Omega)
    %

    %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    % update omega

    somega = somega + dt * ( visc_a + visc_b  + lamb_a );
    
    %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    % check angular momentum conservation


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
