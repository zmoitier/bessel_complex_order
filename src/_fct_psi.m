## Author: Zoïs Moitier

## usage: [psi] = _fct_psi(w, zeta)
##
## Return psi = 2 * ((1-w^2) / (4*zeta) ) ^ (1/4) / w.
##

function [psi] = _fct_psi(w, zeta)
  close_0 = abs(zeta) < 1e-5;
  far_0 = ~close_0;
  psi = zeros(size(zeta));
  
  # zeta close to 0
  psi(close_0) = polyval(
    [5.142857142857142860e-1, 1.007936839915898530e+0, 1.587401051968199470e+0],
    zeta(close_0)
  );
  
  # zeta far form 0
  psi(far_0) = 2 .* power((1 .- w(far_0) .^ 2) ./ (4 .* zeta(far_0)), 0.25) ./ w(far_0);
endfunction
