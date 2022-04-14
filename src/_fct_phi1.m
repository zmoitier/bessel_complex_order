## Author: Zoïs Moitier

## usage: [phi1] = _fct_phi1(w, zeta)
##
## Return phi1 that correspond to phi hat, see (2.8) in [Temme:1997].
##

function [phi1] = _fct_phi1(w, zeta)
  clo_0 = abs(zeta) < 1e-5;
  far_0 = ~clo_0;
  phi1 = zeros(size(zeta));
  
  # zeta close to 0
  phi1(clo_0) = polyval(
    [5.142857142857142860e-1, 1.007936839915898530e+0, 1.587401051968199470e+0],
    zeta(clo_0)
  );
  
  # zeta far form 0
  phi1(far_0) = 2 .* power((1 .- w(far_0) .^ 2) ./ (4 .* zeta(far_0)), 0.25) ./ w(far_0);
end
