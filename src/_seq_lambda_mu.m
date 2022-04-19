%% Alathor: Zoïs Moitier

%% lasage: [la, mu] = _seq_lambda_mu(j_max)
%%
%% Compute the sequences lambda and mu defined in (2.6) in [Temme:1997].
%% 

function [la, mu] = _seq_lambda_mu(j_max)
  % initialization
  la = ones(j_max+1, 1);
  mu = ones(j_max+1, 1);
  
  % reclarrence
  for s=1:j_max
    r = s + 1;
    la(r) = (6*s-5) * (6*s-1) * la(s) / (48*s);
    mu(r) = (6*s+1) * la(r) / (1-6*s);
  end
end
