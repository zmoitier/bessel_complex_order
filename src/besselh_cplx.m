%% Author: Zoïs Moitier

%% usage: [Hvz] = besselh_cplx(nu, k, z)
%%
%% Compute the bessel function H^(k) based on uniform asymptotic expansions for large
%% order describe using the relation https://dlmf.nist.gov/9.2.E11 and in (2.1) in
%% [Temme:1997].
%% 
%% [Temme:1997]
%% N. Temme, Numerical algorithms for uniform Airy-type asymptotic expansions.
%% Numerical Algorithms 15, 207–225 (1997).
%% https://doi.org/10.1023/A:1019197921337
%%

function [Hvz] = besselh_cplx(nu, k, z)
  if ~isscalar(nu) && ~isscalar(z) && ~isequal(size(nu), size(z))
    error("If none of the inputs are scalar, then they must have the same size.")
  end
  
  w = z ./ nu;
  zeta = _fct_zeta(w);
  
  switch k
    case 1
      a0 = 2 * exp(-1i * pi / 3);
      b1 = 2 * exp(1i * pi / 3);
      nu23z = power(nu, 2/3) .* zeta .* exp(2i * pi / 3);
    case 2
      a0 = 2 * exp(1i * pi / 3);
      b1 = 2 * exp(-1i * pi / 3);
      nu23z = power(nu, 2/3) .* zeta .* exp(-2i * pi / 3);
    otherwise
      error("k must be 1 or 2")
  end
  
  [A, B] = _fct_AB(nu, w, zeta, 3);
  Hvz = _fct_phi0(w, zeta) .* (...
    b1 .* airy(1, nu23z) .* B .* power(nu, -5/3) ...
    .+ a0 .* airy(0, nu23z) .* A .* power(nu, -1/3)...
  );
end
