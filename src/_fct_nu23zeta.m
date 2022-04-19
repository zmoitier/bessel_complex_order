%% Author: Zoïs Moitier

%% usage: [nu23z] = _fct_nu23zeta(nu, zeta)
%%
%% Compute the nu ^ (2/3) \* zeta and switch the sign of x - 0i to x + 0i. This is to
%% avoid the Airy bug in Octave, see .
%%

function [nu23z] = _fct_nu23zeta(nu, zeta)
  nu23z = power(nu, 2/3) .* zeta;
  tmp = abs(imag(nu23z)) < (16*eps);
  nu23z(tmp) = (nu23z(tmp) .+ 1i) .- 1i;
end
