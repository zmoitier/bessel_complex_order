## Author: Zoïs Moitier

## usage: [phi] = _fct_phi(w, zeta)
##
## Return phi = (4*zeta / (1-w^2)) ^ (1/4).
##

function [phi] = _fct_phi(w, zeta)
  close_0 = abs(zeta) < 1e-6;
  far_0 = ~close_0;
  phi = zeros(size(zeta));
  
  # zeta close to 0
  phi(close_0) = polyval(
    [2.040944209673399320e-2, 2.000000000000000000e-1, 1.259921049894873160e+0],
    zeta(close_0)
  );
  
  # zeta far form 0
  phi(far_0) = power(4 .* zeta(far_0) ./ (1 .- w(far_0) .^ 2), 0.25);
endfunction
