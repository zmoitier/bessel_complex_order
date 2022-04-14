clear;
addpath("src/");
hold on;

N = 128;
x = logspace(-16, 0, N);

t = 1; #exp(1i*pi/3)
w0 = 1 .- x .* t;
w2 = 1 .+ x .* t;

z0 = abs(_fct_zeta(w0));
z2 = abs(_fct_zeta(w2));

loglog(x, z0)
loglog(x, z2)

##N = 64;
##x = linspace(-2e-8, 2e-8, N);
##W = 1 .+ x + 1i .* x';
##Z = _fct_zeta(W);
##
##subplot(1, 2, 1);
##im = pcolor(x, x, real(Z));
##set(im, 'EdgeColor', 'none');
##colorbar();
##
##subplot(1, 2, 2);
##im = pcolor(x, x, imag(Z));
##set(im, 'EdgeColor', 'none');
##colorbar();

hold off;
