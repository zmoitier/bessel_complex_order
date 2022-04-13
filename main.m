clear;
addpath("src/");
hold on;

N = 128;
x = linspace(0, 2, N);
y = linspace(-1, 1, N);
[X, Y] = meshgrid(x, y);
W = X .+ 1i .* Y;

Z = _fct_zeta(W);

zs1 = abs(Z) <= 1;
zg1 = ~zs1;

F = zeros(size(Z));
F = sqrt(Z);
#F(zg1) = (-1i) .* sqrt(-Z(zg1));

subplot(1, 2, 1);
im = pcolor(X, Y, real(F));
set(im, 'EdgeColor', 'none');

subplot(1, 2, 2);
im = pcolor(X, Y, imag(F));
set(im, 'EdgeColor', 'none');

hold off;
