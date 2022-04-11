## Author: Zoïs Moitier

## usage: [u, v] = _seq_uv(j_max)
##
## Compute the sequence (3/2)^j*u_j and (3/2)^j*v_j where u_j and v_j are defined in [NIST:9.7.2].
## 

function [u, v] = _seq_uv(j_max)
  % initialization
  u = ones(j_max+1, 1);
  v = ones(j_max+1, 1);
  
  % recurrence
  for s=1:j_max
    r = s + 1;
    u(r) = (6*s-5) * (6*s-1) * u(s) / (48*s);
    v(r) = (6*s+1) * u(r) / (1-6*s);
  end
end
