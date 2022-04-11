## Author: Zoïs Moitier

## usage: [phi] = _fct_phi(w, zeta)
##
## Get the value of phi from w and zeta using (2.4) in [NIST:10.20(i)].
##

function [phi] = _fct_phi(w, zeta)
  phi = power(4 .* zeta ./ (1 .- w .^ 2), 0.25);
  
  # w = 1
  is_one = abs(w - 1) < 1e-14;
  phi(is_one) = power(2, 1/3);
endfunction
