function omega = om_b_trans(somega)
  global jmax leg_n1_si
  for j=1:jmax
    omega(:,j) = leg_n1_si' * somega(:,j);
  end
