## Author: Zoïs Moitier

## usage: [psi] = _fct_psi(w, zeta)
##
## Return psi = 2 * ((1-w^2) / (4*zeta) ) ^ (1/4) / w.
##

function [psi] = _fct_psi(w, zeta)
  psi = 2 .* power((1 .- w .^ 2) ./ (4 .* zeta), 0.25) ./ w;
  
  # w = 1
  is_one = abs(w - 1) < 1e-14;
  psi(is_one) = power(2, 2/3);
endfunction
