%% Author: Zoïs Moitier

%% usage: [U] = _poly_U(j_max)
%%
%% Compute the polynoms U from 0 to j_max from (2.7) in [Temme:1997].
%%

function U = _poly_U(j_max)
  D = [-0.5, 0, 0.5, 0, 0];
  I = [-0.625, 0, 0.125];
  
  % initialization
  U{1} = [1];
  
  % recurrence
  if j_max >= 1
    U{2} = polyint(I);
  end
  
  for j=3:(j_max+1)
    U{j} = conv(D, polyder(U{j-1})) + polyint(conv(I, U{j-1}));
  end
end
