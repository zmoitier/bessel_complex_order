%% Author: Zoïs Moitier

%% usage: [zeta] = _fct_zeta(w)
%%
%% Get the value zeta from w using (2.3) in [Temme:1997].
%%

function [zeta] = _fct_zeta(w)
  tp = abs(w .- 1) < 1e-8;
  inn = (abs(w) <= 1) & ~tp;
  out = ~inn & ~tp;
  
  zeta = zeros(size(w));
  
  % |w| <= 1
  q = sqrt(1 .- w(inn) .^ 2);
  zeta(inn) = power(1.5 .* (log((1 .+ q) ./ w(inn)) .- q), 2/3);
  
  % turning point
  zeta(tp) = polyval([3.779763149684619490e-1, -1.259921049894873160e+0, 0], w(tp) .- 1);
  
  % |w| > 1
  zeta(out) = -power(1.5 .* (sqrt(w(out) .^ 2 .- 1) .- asec(w(out))), 2/3);
end
