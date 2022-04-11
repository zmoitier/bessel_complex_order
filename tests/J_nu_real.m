clear;
addpath("../src");


N = 124;
nu = logspace(0, 2, N);
w = linspace(0.5, 1.5, N);

[NU, W] = meshgrid(nu, w);
Z = NU .* W;

Jvz = besselc_J(NU, Z);
ref = besselj(NU, Z);
err = log10(abs(Jvz ./ ref .- 1));

## 

hold on;
h = pcolor(log10(NU), W, err);
## shading interp;
set(h, 'EdgeColor', 'none');
axis tight;
colorbar();

plot([log10(nu(1)), log10(nu(N))], [0.98, 0.98], "k")
plot([log10(nu(1)), log10(nu(N))], [1.02, 1.02], "k")

## p = polyfit(log10(z), log10(err), 1);
## disp(["slope = ", num2str(p(1))])

hold off;
