%% Author: Zoïs Moitier

%% usage: [A, B] = _fct_AB_far_tp(nu, w, zeta, k_max)
%%
%% Get the value for the functions A_nu(zeta) and B_nu(zeta) see (2.2) and (2.5) in
%% [Temme:1997]. This is numericaly stable if zeta is not to close of 0.
%%

function [A, B] = _fct_AB_far_tp(nu, w, zeta, k_max)
  [la, mu] = _seq_lambda_mu(2*k_max+1);
  U = _poly_U(2*k_max+1);
  
  ws1 = abs(w) <= 1;
  wg1 = abs(w) > 1;
  
  p = zeros(size(w));
  p(ws1) = power(1 .- w(ws1) .^ 2, -0.5);
  p(wg1) = 1i .* power(w(wg1) .^ 2 .- 1, -0.5);
  
  zeta32 = zeros(size(zeta));
  zeta32(ws1) = power(zeta(ws1), -1.5);
  zeta32(wg1) = (-1i) .* power(-zeta(wg1), -1.5);
  
  zeta12 = zeros(size(zeta));
  zeta12(ws1) = power(zeta(ws1), -0.5);
  zeta12(wg1) = 1i .* power(-zeta(wg1), -0.5);
  
  inv_nu2 = power(nu, -2);
  
  A = _coef_a(p, zeta32, mu, U, k_max);
  B = _coef_b(p, zeta32, la, U, k_max);
  for k = (k_max-1):(-1):0
    A = _coef_a(p, zeta32, mu, U, k) .+ inv_nu2 .* A;
    B = _coef_b(p, zeta32, la, U, k) .+ inv_nu2 .* B;
  end
  B .*= -zeta12;
end

function [Ak] = _coef_a(p, zeta32, mu, U, k)
  Ak = polyval(mu(2*k+1) .* U{1}, p);
  for j = (2*k-1):(-1):0
    Ak = polyval(mu(j+1) .* U{2*k-j+1}, p) .+ zeta32 .* Ak;
  end
end

function [Bk] = _coef_b(p, zeta32, la, U, k)
  Bk = polyval(la(2*k+2) .* U{1}, p);
  for j = (2*k):(-1):0
    Bk = polyval(la(j+1) .* U{2*k-j+2}, p) .+ zeta32 .* Bk;
  end
end
