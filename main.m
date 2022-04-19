clear;
addpath("src");
## hold on;

N = 3;
nu = 2 .* logspace(0, 2, N);

r = 0.5 .* linspace(-1, 1, N);
w = 1 .+ r .* exp(1i * pi / 6);

[NU, W] = meshgrid(nu, w);
Z = NU .* W;

##disp(num2str(NU(62, end), "%.18f"))
##disp(num2str(Z(62, end), "%.18f"))
##
##J0 = (besselc_J(NU .- 1, Z) .- besselc_J(NU .+ 1, Z)) ./ 2;
##J1 = besselc_Jp(NU, Z);
##
##disp(["J0 = ", num2str(J0(62, end), "%.18e")])
##disp(["J1 = ", num2str(J1(62, end), "%.18e")])

disp("")

##disp("Matrix version")
besselc_J(NU .- 1, Z)(2, end);
##disp(besselc_J(NU .- 1, Z)(62, end))
##disp(besselj(NU .- 1, Z)(62, end))

##disp("Scalar version")
besselc_J(NU(2, end) .- 1, Z(2, end));
##disp(besselc_J(NU(62, end) .- 1, Z(62, end)))
##disp(besselj(NU(62, end) .- 1, Z(62, end)))

disp("")

## hold off;


airy(0, z, true) - airy(0, z) * exp(2/3 * z * sqrt(z))
airy(1, z, true) - airy(1, z) * exp(2/3 * z * sqrt(z))
airy(2, z, true) - airy(2, z) * exp(-abs(real(2/3 * z * sqrt(z))))
airy(3, z, true) - airy(3, z) * exp(-abs(real(2/3 * z * sqrt(z))))