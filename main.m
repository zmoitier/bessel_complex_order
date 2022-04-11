clear;
addpath("src/");
hold on;

z = 1.5 .* exp(1i * linspace(0, 2*pi, 129)(1:end-1));
zeta = _fct_zeta(z);

cond = z < 1;
z(cond)

e1 = power(1 .- z .^ 2, -0.5);
e2 = 1i .* power(z .^ 2 .- 1, -0.5);

plot(real(e1), imag(e1), "+");
plot(real(e2), imag(e2), "x");

axis("equal")

hold off;
