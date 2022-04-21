%% Author: Zoïs Moitier

%% usage: [Jvz] = besselj_cplx(nu, z)
%%
%% Compute the bessel function J based on uniform asymptotic expansions for large order
%% describe in (2.1) in [Temme:1997].
%% 
%% [Temme:1997]
%% N. Temme, Numerical algorithms for uniform Airy-type asymptotic expansions.
%% Numerical Algorithms 15, 207–225 (1997).
%% https://doi.org/10.1023/A:1019197921337
%%

function [Jvz] = besselj_cplx(nu, z)
  if ~isscalar(nu) && ~isscalar(z) && ~isequal(size(nu), size(z))
    error("If none of the inputs are scalar, then they must have the same size.")
  end
  
  w = z ./ nu;
  zeta = _fct_zeta(w);
  
  nu23z = _fct_nu23zeta(nu, zeta);
  [A, B] = _fct_AB(nu, w, zeta, 3);
  Jvz = _fct_phi0(w, zeta) .* (...
    airy(1, nu23z) .* B .* power(nu, -5/3) .+ airy(0, nu23z) .* A .* power(nu, -1/3)...
  );
end
