%% Author: Zoïs Moitier

%% usage: [A, B] = _fct_AB(nu, w, zeta, k_max)
%%
%% Get the value for the function A_nu(zeta) and B_nu(zeta).
%%

function [A, B] = _fct_AB(nu, w, zeta, k_max)
  [clo_tp, far_tp] = _index_tp(zeta, 0.175);
  
  A = zeros(size(zeta));
  B = zeros(size(zeta));
  if isscalar(nu)
    if ~isempty(clo_tp)
      [A(clo_tp), B(clo_tp)] = _fct_AB_clo_tp(nu, zeta(clo_tp), 3);
    end
    if ~isempty(far_tp)
      [A(far_tp), B(far_tp)] = _fct_AB_far_tp(nu, w(far_tp), zeta(far_tp), 3);
    end
  else
    if ~isempty(clo_tp)
      [A(clo_tp), B(clo_tp)] = _fct_AB_clo_tp(nu(clo_tp), zeta(clo_tp), 3);
    end
    if ~isempty(far_tp)
      [A(far_tp), B(far_tp)] = _fct_AB_far_tp(nu(far_tp), w(far_tp), zeta(far_tp), 3);
    end
  end
end
