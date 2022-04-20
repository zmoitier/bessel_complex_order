%% Author: Zoïs Moitier

%% usage: [Jpvz] = besselc_Jp(nu, z)
%%
%% Compute the derivative of bessel function J' based on uniform asymptotic expansions
%% for large order describe in (2.8) in [Temme:1997].
%% 
%% [Temme:1997]
%% N. Temme, Numerical algorithms for uniform Airy-type asymptotic expansions.
%% Numerical Algorithms 15, 207–225 (1997).
%% https://doi.org/10.1023/A:1019197921337
%%

function [Jpvz] = besselc_Jp(nu, z)
  w = z ./ nu;
  zeta = _fct_zeta(w);
  
  clo_tp = abs(zeta) < 0.175;
  far_tp = ~clo_tp;
  
  C = zeros(size(zeta));
  D = zeros(size(zeta));
  [C(clo_tp), D(clo_tp)] = _fct_CD_tp(nu(clo_tp), zeta(clo_tp), 3);
  [C(far_tp), D(far_tp)] = _fct_CD(nu(far_tp), w(far_tp), zeta(far_tp), 3);
  
  nu23z = _fct_nu23zeta(nu, zeta);
  Jpvz = -_fct_phi1(w, zeta) .* (
    airy(0, nu23z) .* C .* power(nu, -4/3) .+ airy(1, nu23z) .* D .* power(nu, -2/3)
  );
end
