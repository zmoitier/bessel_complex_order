clear;
addpath("../src");

N = 124;
nu = logspace(0, 2, N);
w = linspace(0.5, 1.5, N);

[NU, W] = meshgrid(nu, w);
Z = NU .* W;

Jvz = besselc_J(NU, Z);
ref = besselj(NU, Z);
err = abs(Jvz ./ ref .- 1);
err(err <= 0) = 1e-16;

wmin = 0.95;
wmax = 1.05;

############
#### 2D plot
figure(1);
hold on;

h = pcolor(log10(NU), W, log10(err));
set(h, 'EdgeColor', 'none');

xlabel("log_{10}(\\nu)", "fontsize", 16)
ylabel("w", "fontsize", 16)
axis tight;
colorbar();

title("log_{10} of the relative error for J_\\nu(\\nu w)", "fontsize", 16)

hold off;

#####################
#### nu -> J_nu(nu w)
figure(2);
hold on;

w_vec = [0.9*wmin, wmin, 1, wmax, 1.1*wmax];
leg_name = {};
for j = 1:length(w_vec)
  [_ i] = min(abs(w .- w_vec(j)));
  loglog(nu, err(i, :), "+--");
  leg_name{j} = ["w = ", num2str(w(i))];
endfor

xlabel("\\nu", "fontsize", 16);
ylabel("relative error", "fontsize", 16);
axis tight;
grid on;
legend(leg_name, "fontsize", 16);

title("log_{10} of the relative error for \\nu -> J_\\nu(\\nu w)", "fontsize", 16)

hold off;

####################
#### w -> J_nu(nu w)
figure(3);
hold on;

[_ i0] = min(abs(nu .- 1));
loglog(w, err(:, i0), "+--");

[_ i1] = min(abs(nu .- 10));
loglog(w, err(:, i1), "+--");

[_ i2] = min(abs(nu .- 100));
loglog(w, err(:, i2), "+--");

xlabel("w", "fontsize", 16)
ylabel("relative error", "fontsize", 16)
axis tight;
grid on;
legend(
  {["\\nu = ", num2str(nu(i0))], ["\\nu = ", num2str(nu(i1))], ["\\nu = ", num2str(nu(i2))]},
  "fontsize", 16
)

title("log_{10} of the relative error for w -> J_\\nu(\\nu w)", "fontsize", 16)

hold off;