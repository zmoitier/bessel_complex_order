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
  
  [A, B] = _fct_AB(nu, w, zeta, 2);
  
  nu23z = power(nu, 2/3) .* zeta;
  Jvz = phi .* (
      airy(0, nu23z) .* A .* power(nu, -1/3) .+ airy(1, nu23z) .* B .* power(nu, -5/3)
    );
endfunction
