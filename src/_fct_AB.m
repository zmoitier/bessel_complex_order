## Author: Zoïs Moitier

## usage: [A, B] = _fct_AB(nu, w, zeta)
##
## Get the value for the function A_nu(zeta) and B_nu(zeta) used in [NIST:10.20.4-6].
##

function [A, B] = _fct_AB(nu, w, zeta, k_max)
  [u, v] = _seq_uv(2*k_max+1);
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
  
  A = _coef_A(p, zeta32, v, U, k_max);
  B = _coef_B(p, zeta32, u, U, k_max);
  for k = (k_max-1):(-1):0
    A = _coef_A(p, zeta32, v, U, k) .+ inv_nu2 .* A;
    B = _coef_B(p, zeta32, u, U, k) .+ inv_nu2 .* B;
  endfor
  B .*= -zeta12;
endfunction

function [Ak] = _coef_A(p, zeta32, v, U, k)
  Ak = polyval(v(2*k+1) .* U{1}, p);
  for j = (2*k-1):(-1):0
    Ak = polyval(v(j+1) .* U{2*k-j+1}, p) .+ zeta32 .* Ak;
  endfor
endfunction

function [Bk] = _coef_B(p, zeta32, u, U, k)
  Bk = polyval(u(2*k+2) .* U{1}, p);
  for j = (2*k):(-1):0
    Bk = polyval(u(j+1) .* U{2*k-j+2}, p) .+ zeta32 .* Bk;
  endfor
endfunction
