%% Author: Zoïs Moitier

%% usage: [Ypvz] = besselyp_cplx(nu, z)
%%
%% Compute the derivative of bessel function J' based on uniform asymptotic expansions
%% for large order describe in (2.8) in [Temme:1997].
%% 
%% [Temme:1997]
%% N. Temme, Numerical algorithms for uniform Airy-type asymptotic expansions.
%% Numerical Algorithms 15, 207–225 (1997).
%% https://doi.org/10.1023/A:1019197921337
%%

function [Ypvz] = besselyp_cplx(nu, z)
  if ~isscalar(nu) && ~isscalar(z) && ~isequal(size(nu), size(z))
    error("If none of the inputs are scalar, then they must have the same size.")
  end
  
  w = z ./ nu;
  zeta = _fct_zeta(w);
  
  nu23z = power(nu, 2/3) .* zeta;
  [C, D] = _fct_CD(nu, w, zeta, 3);
  Ypvz = _fct_phi1(w, zeta) .* (...
    airy(2, nu23z) .* C .* power(nu, -4/3) .+ airy(3, nu23z) .* D .* power(nu, -2/3)...
  );
end
