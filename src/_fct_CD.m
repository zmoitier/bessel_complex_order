%% Author: Zoïs Moitier

%% usage: [C, D] = _fct_CD(nu, w, zeta, k_max)
%%
%% Get the value for the function C_nu(zeta) and D_nu(zeta).
%%

function [C, D] = _fct_CD(nu, w, zeta, k_max)
  [clo_tp, far_tp] = _index_tp(zeta, 0.175);
  
  C = zeros(size(zeta));
  D = zeros(size(zeta));
  if isscalar(nu)
    if ~isempty(clo_tp)
      [C(clo_tp), D(clo_tp)] = _fct_CD_clo_tp(nu, zeta(clo_tp), 3);
    end
    if ~isempty(far_tp)
      [C(far_tp), D(far_tp)] = _fct_CD_far_tp(nu, w(far_tp), zeta(far_tp), 3);
    end
  else
    if ~isempty(clo_tp)
      [C(clo_tp), D(clo_tp)] = _fct_CD_clo_tp(nu(clo_tp), zeta(clo_tp), 3);
    end
    if ~isempty(far_tp)
      [C(far_tp), D(far_tp)] = _fct_CD_far_tp(nu(far_tp), w(far_tp), zeta(far_tp), 3);
    end
  end
end
