function somega = om_f_trans(omega)
% This is the lousy cheapo way of doing the forward transform,
% which really ought to be done using Gauss-Legendre, but I haven't
% bothered since it doesn't really affect the core algorithm.
  global nmax imax norm leg_n1 s1d
  dth = pi/(imax-1),  % lo
  for n=1:nmax,
    %   for j=1:jmax,  % loop over r
    %   Note that P_1^1(cos(th)) / sin(th) = -1, so this takes care of constant.
    somega(n,:) = norm(n) * (leg_n1(n,:) .* s1d.^2) * omega  * dth;
    %      somega(n,j) = norm(n) * (leg_n1(n,:) .* s1d.^2) * omega(:,j);
    %   end			   
  end
