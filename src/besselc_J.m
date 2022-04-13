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
  
  clo_tp = abs(zeta) < 0.1;
  far_tp = ~clo_tp;
  
  A = zeros(size(zeta));
  B = zeros(size(zeta));
  [A(clo_tp), B(clo_tp)] = _fct_AB_tp(nu(clo_tp), zeta(clo_tp));
  [A(far_tp), B(far_tp)] = _fct_AB(nu(far_tp), w(far_tp), zeta(far_tp), 2);
  
  nu23z = power(nu, 2/3) .* zeta;
  Jvz = phi .* (
    airy(0, nu23z) .* A .* power(nu, -1/3) .+ airy(1, nu23z) .* B .* power(nu, -5/3)
  );
endfunction
