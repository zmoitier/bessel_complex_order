## Author: Zoïs Moitier

## usage: [Jpvz] = besselc_Jp(nu, z)
##
## Compute the bessel function J' based on uniform asymptotic expansions for large order
## describe in [NIST:10.20].
##

function [Jpvz] = besselc_Jp(nu, z)
  w = z ./ nu;
  zeta = _fct_zeta(w);
  psi = _fct_psi(w, zeta);
  
  [C, D] = _fct_CD(nu, w, zeta, 2);
  
  nu23z = power(nu, 2/3) .* zeta;
  Jpvz = -psi .* (
      airy(0, nu23z) .* C .* power(nu, -4/3) .+ airy(1, nu23z) .* D .* power(nu, -2/3)
    );
endfunction
