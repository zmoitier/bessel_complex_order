%% Author: Zoïs Moitier

%% usage: [Yvz] = bessely_cplx(nu, z)
%%
%% Compute the bessel function Y based on uniform asymptotic expansions for large order
%% describe in (2.1) in [Temme:1997].
%% 
%% [Temme:1997]
%% N. Temme, Numerical algorithms for uniform Airy-type asymptotic expansions.
%% Numerical Algorithms 15, 207–225 (1997).
%% https://doi.org/10.1023/A:1019197921337
%%

function [Yvz] = bessely_cplx(nu, z)
  if ~isscalar(nu) && ~isscalar(z) && ~isequal(size(nu), size(z))
    error("If none of the inputs are scalar, then they must have the same size.")
  end
  
  w = z ./ nu;
  zeta = _fct_zeta(w);
  
  nu23z = power(nu, 2/3) .* zeta;
  [A, B] = _fct_AB(nu, w, zeta, 3);
  Yvz = -_fct_phi0(w, zeta) .* (
    airy(3, nu23z) .* B .* power(nu, -5/3) .+ airy(2, nu23z) .* A .* power(nu, -1/3)
  );
end
