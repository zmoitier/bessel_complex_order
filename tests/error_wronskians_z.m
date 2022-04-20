clear;
addpath("../src");

nu = 100 * exp(-1i * pi / 6);

N = 127;
x = real(nu) .+ linspace(-25, 25, N);
y = imag(nu) .+ linspace(-25, 25, N);
[X, Y] = meshgrid(x, y);
Z = X .+ 1i .* Y;

NU = nu .* ones(N, N);
J0 = besselc_J(NU, Z);
J1 = besselc_Jp(NU, Z);
Y0 = besselc_Y(NU, Z);
Y1 = besselc_Yp(NU, Z);

Wron = (Y1 ./ Y0 .- J1 ./ J0) .* J0 .* Y0;

err = abs(Wron .* Z .* (pi/2) .- 1);
%err(err <= 1e-16) = 1e-16;

figure(1);
hold on;

h = pcolor(x, y, log10(err));
set(h, 'EdgeColor', 'none');

xlabel("\\Re(z)", "fontsize", 16)
ylabel("\\Im(z)", "fontsize", 16)
axis tight;
colorbar();
%caxis ([-16 0]);

title("log_{10} of the relative error for W(J_\\nu, Y_\\nu)(z)", "fontsize", 16)

hold off;