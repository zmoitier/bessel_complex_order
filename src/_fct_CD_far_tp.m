%% Author: Zoïs Moitier

%% usage: [C, D] = _fct_CD_far_tp(nu, w, zeta, k_max)
%%
%% Get the value for the functions C_nu(zeta) and D_nu(zeta) see (2.10) and (2.12) in
%% [Temme:1997]. This is numericaly stable if zeta is not to close of 0.
%%

function [C, D] = _fct_CD_far_tp(nu, w, zeta, k_max)
  [la, mu] = _seq_lambda_mu(2*k_max+1);
  V = _poly_V(2*k_max+1);
  
  ws1 = abs(w) <= 1;
  wg1 = abs(w) > 1;
  
  p = zeros(size(w));
  p(ws1) = power(1 .- w(ws1) .^ 2, -0.5);
  p(wg1) = 1i .* power(w(wg1) .^ 2 .- 1, -0.5);
  
  zeta32 = zeros(size(zeta));
  zeta32(ws1) = power(zeta(ws1), -1.5);
  zeta32(wg1) = (-1i) .* power(-zeta(wg1), -1.5);
  
  zeta12 = zeros(size(zeta));
  zeta12(ws1) = sqrt(zeta(ws1));
  zeta12(wg1) = (-1i) .* sqrt(-zeta(wg1));
  
  inv_nu2 = power(nu, -2);
  
  C = _coef_C(p, zeta32, mu, V, k_max);
  D = _coef_D(p, zeta32, la, V, k_max);
  for k = (k_max-1):(-1):0
    C = _coef_C(p, zeta32, mu, V, k) .+ inv_nu2 .* C;
    D = _coef_D(p, zeta32, la, V, k) .+ inv_nu2 .* D;
  end
  C = -zeta12 .* C;
end

function [Ck] = _coef_C(p, zeta32, mu, V, k)
  Ck = polyval(mu(2*k+2) .* V{1}, p);
  for j = (2*k):(-1):0
    Ck = polyval(mu(j+1) .* V{2*k-j+2}, p) .+ zeta32 .* Ck;
  end
end

function [Dk] = _coef_D(p, zeta32, la, V, k)
  Dk = polyval(la(2*k+1) .* V{1}, p);
  for j = (2*k-1):(-1):0
    Dk = polyval(la(j+1) .* V{2*k-j+1}, p) .+ zeta32 .* Dk;
  end
end
