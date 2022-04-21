clear;
addpath("../src");

nu = 400 - 1i;

N = 127;
x = real(nu) .+ linspace(-20, 20, N);
y = imag(nu) .+ linspace(-20, 20, N);
[X, Y] = meshgrid(x, y);
Z = X .+ 1i .* Y;

J0 = besselj_cplx(nu, Z);
J1 = besseljp_cplx(nu, Z);
Y0 = bessely_cplx(nu, Z);
Y1 = besselyp_cplx(nu, Z);
WB = (Y1 ./ Y0 .- J1 ./ J0) .* J0 .* Y0;
eB = abs(WB .* Z .* (pi/2) .- 1);

H10 = besselh_cplx(nu, 1, Z);
H11 = besselhp_cplx(nu, 1, Z);
H20 = besselh_cplx(nu, 2, Z);
H21 = besselhp_cplx(nu, 2, Z);
WH = (H21 ./ H20 .- H11 ./ H10) .* H20 .* H10;
eH = abs(WH .* Z .* (-pi/(4i)) .- 1);

figure(1);

subplot(1, 2, 1);
h = pcolor(x, y, log10(eB));
set(h, 'EdgeColor', 'none');

xlabel("\\Re(z)", "fontsize", 16)
ylabel("\\Im(z)", "fontsize", 16)
axis tight;
colorbar();
caxis ([-16 0]);

title("log_{10} of the relative error for W(J_\\nu, Y_\\nu)(z)", "fontsize", 16)

subplot(1, 2, 2);
h = pcolor(x, y, log10(eH));
set(h, 'EdgeColor', 'none');

xlabel("\\Re(z)", "fontsize", 16)
ylabel("\\Im(z)", "fontsize", 16)
axis tight;
colorbar();
caxis ([-16 0]);

title("log_{10} of the relative error for W(H^{(1)}_\\nu, H^{(2)}_\\nu)(z)", "fontsize", 16)
