## Author: Zoïs Moitier

## usage: [V] = _poly_V(j_max)
##
## Compute the polynoms V from 0 to j_max from (2.7) in [NIST:10.41.9].
##

function V = _poly_V(j_max)
  U = _poly_U(j_max)
  D0 = [-0.5, 0, 0.5, 0];
  D1 = [-1, 0, 1, 0, 0];
  
  # initialization
  V{1} = [1];
  
  # recurrence
  if j_max >= 1
    V{2} = U{2} - D0;
  endif
  
  for j=2:(j_max+1)
    V{j+1} = U{j+1} - conv(D0, U{j}) - conv(D1, polyder(U{j}));
  endfor
endfunction
