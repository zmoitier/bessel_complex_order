## Author: Zoïs Moitier

## usage: [phi0] = _fct_phi0(w, zeta)
##
## Return phi0 = (4*zeta / (1-w^2)) ^ (1/4), see (2.4) in [Temme:1997].
##

function [phi0] = _fct_phi0(w, zeta)
  clo_0 = abs(zeta) < 1e-6;
  far_0 = ~clo_0;
  phi0 = zeros(size(zeta));
  
  # zeta close to 0
  phi0(clo_0) = polyval(
    [2.040944209673399320e-2, 2.000000000000000000e-1, 1.259921049894873160e+0],
    zeta(clo_0)
  );
  
  # zeta far form 0
  phi0(far_0) = power(4 .* zeta(far_0) ./ (1 .- w(far_0) .^ 2), 0.25);
end
