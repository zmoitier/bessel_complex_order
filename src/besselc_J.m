## Author: Zoïs Moitier

## usage: [Jvz] = besselc_J(nu, z)
##
## Compute the bessel function J based on uniform asymptotic expansions for large order
## describe in (2.1) in [Temme:1997].
## 
## [Temme:1997]
## N. Temme, Numerical algorithms for uniform Airy-type asymptotic expansions.
## Numerical Algorithms 15, 207–225 (1997).
## https://doi.org/10.1023/A:1019197921337
##

function [Jvz] = besselc_J(nu, z)
  w = z ./ nu;
  zeta = _fct_zeta(w);
  
  clo_tp = abs(zeta) < 0.175;
  far_tp = ~clo_tp;
  
  A = zeros(size(zeta));
  B = zeros(size(zeta));
  [A(clo_tp), B(clo_tp)] = _fct_AB_tp(nu(clo_tp), zeta(clo_tp), 3);
  [A(far_tp), B(far_tp)] = _fct_AB(nu(far_tp), w(far_tp), zeta(far_tp), 3);
  
  nu23z = power(nu, 2/3) .* zeta;
  Jvz = _fct_phi0(w, zeta) .* (
    airy(0, nu23z) .* A .* power(nu, -1/3) .+ airy(1, nu23z) .* B .* power(nu, -5/3)
  );
end
