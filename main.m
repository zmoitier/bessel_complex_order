clear;
addpath("src/");
hold on;

w = horzcat(linspace(0.1, 1, 64), linspace(1, 1.9, 64));
zeta = _fct_zeta(w);
phi = _fct_phi(w, zeta);
psi = _fct_psi(w, zeta);

plot(w, phi);
plot(w, psi);

hold off;
