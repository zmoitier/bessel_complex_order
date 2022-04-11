## Author: Zoïs Moitier

## usage: [zeta] = _fct_zeta(w)
##
## Get the value zeta from w using (2.3) in [NIST:10.20(i)].
##

function [zeta] = _fct_zeta(w)
  zeta = zeros(size(w));
  
  # |w| <= 1
  is_smaller = abs(w) <= 1;
  q = sqrt(1 .- w(is_smaller) .^ 2);
  zeta(is_smaller) = power(1.5 .* (log((1 .+ q) ./ w(is_smaller)) .- q), 2/3);
  
  # |w| > 1
  is_larger = abs(w) > 1;
  zeta(is_larger) = -power(1.5 .* (sqrt(w(is_larger) .^ 2 .- 1) .- asec(w(is_larger))), 2/3);
endfunction
