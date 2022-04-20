clear;
addpath("../src");

z = 100;

N = 127;
x = real(z) .+ linspace(-25, 25, N);
y = imag(z) .+ linspace(-25, 25, N);
[X, Y] = meshgrid(x, y);
NU = X .+ 1i .* Y;

J0 = besselj_cplx(NU, z);
J1 = besseljp_cplx(NU, z);
Y0 = bessely_cplx(NU, z);
Y1 = besselyp_cplx(NU, z);
WB = (Y1 ./ Y0 .- J1 ./ J0) .* J0 .* Y0;
eB = abs(WB .* z .* (pi/2) .- 1);

H10 = besselh_cplx(NU, 1, z);
H11 = besselhp_cplx(NU, 1, z);
H20 = besselh_cplx(NU, 2, z);
H21 = besselhp_cplx(NU, 2, z);
WH = (H21 ./ H20 .- H11 ./ H10) .* H20 .* H10;
eH = abs(WH .* z .* (-pi/(4i)) .- 1);

figure(1);

subplot(1, 2, 1);
h = pcolor(x, y, log10(eB));
set(h, 'EdgeColor', 'none');

xlabel("\\Re(\\nu)", "fontsize", 16)
ylabel("\\Im(\\nu)", "fontsize", 16)
axis tight;
colorbar();
caxis ([-16 0]);

title("log_{10} of the relative error for W(J_\\nu, Y_\\nu)(z)", "fontsize", 16)

subplot(1, 2, 2);
h = pcolor(x, y, log10(eH));
set(h, 'EdgeColor', 'none');

xlabel("\\Re(\\nu)", "fontsize", 16)
ylabel("\\Im(\\nu)", "fontsize", 16)
axis tight;
colorbar();
caxis ([-16 0]);

title("log_{10} of the relative error for W(H^{(1)}_\\nu, H^{(2)}_\\nu)(z)", "fontsize", 16)
