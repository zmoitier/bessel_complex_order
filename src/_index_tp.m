%% Author: Zoïs Moitier

%% usage: [clo_tp, far_tp] = _index_tp(zeta, cutoff)
%%
%% Return the index for which the element of zeta are smaller and bigger in modulus to a
%% cutoff.
%%

function [clo_tp, far_tp] = _index_tp(zeta, cutoff)
  is_smaller = abs(zeta) <= cutoff;
  clo_tp = find(is_smaller);
  far_tp = find(~is_smaller);
end
