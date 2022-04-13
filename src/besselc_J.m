## Author: Zoïs Moitier

## usage: [Jvz] = besselc_J(nu, z)
##
## Compute the bessel function J based on uniform asymptotic expansions for large order
## describe in [NIST:10.20].
##

function [Jvz] = besselc_J(nu, z)
  w = z ./ nu;
  zeta = _fct_zeta(w);
  phi = _fct_phi(w, zeta);
  
  ##  z0 = abs(nu) < 20;
  ##  nu50 = abs(nu) < 50;
  ##  z1 = ~z0 & nu50;
  ##  z2 = ~nu50;
  ##  
  ##  A = zeros(size(w));
  ##  B = zeros(size(w));
  ##  [A(z0), B(z0)] = _fct_AB(nu(z0), w(z0), zeta(z0), 4);
  ##  [A(z1), B(z1)] = _fct_AB(nu(z1), w(z1), zeta(z1), 3);
  ##  [A(z2), B(z2)] = _fct_AB(nu(z2), w(z2), zeta(z2), 2);
  
  [A, B] = _fct_AB(nu, w, zeta, 5);
  
  nu23z = power(nu, 2/3) .* zeta;
  Jvz = phi .* (
    airy(0, nu23z) .* A .* power(nu, -1/3) .+ airy(1, nu23z) .* B .* power(nu, -5/3)
  );
endfunction
