clear;
addpath("../src");

N = 127;
nu = logspace(0, 2, N);

r = 0.5 .* linspace(-1, 1, N);
w = 1 .+ r .* exp(1i * pi / 6);

[NU, W] = meshgrid(nu, w);
Z = NU .* W;

H10 = besselh_cplx(NU, 1, Z);
H11 = besselhp_cplx(NU, 1, Z);
H20 = besselh_cplx(NU, 2, Z);
H21 = besselhp_cplx(NU, 2, Z);

Wron = (H21 ./ H20 .- H11 ./ H10) .* H20 .* H10;

err = abs(Wron .* Z .* (-pi/(4i)) .- 1);
err(err <= 1e-16) = 1e-16;

%%%%%%%%%%%%
%%%% 2D plot
figure(1);
hold on;

h = pcolor(log10(abs(nu)), r, log10(err));
set(h, 'EdgeColor', 'none');

xlabel("log_{10}(\\nu)", "fontsize", 16)
ylabel("r", "fontsize", 16)
axis tight;
colorbar();
caxis ([-16 0]);

title("log_{10} of the relative error for W(H^{(1)}_\\nu, H^{(2)}_\\nu)(\\nu w)", "fontsize", 16)

hold off;

%%%%%%%%%%%%%%%%%%%%%
%%%% nu -> J_nu(nu w)
figure(2);
hold on;

r_vec = [-0.25, -0.05, 0, 0.05, 0.25];
leg_name = {};
for j = 1:length(r_vec)
  [_ i] = min(abs(r .- r_vec(j)));
  loglog(abs(nu), err(i, :), "+--");
  leg_name{j} = ["r = ", num2str(r(i))];
endfor

xlabel("\\nu", "fontsize", 16);
ylabel("relative error", "fontsize", 16);
axis tight;
grid on;
legend(leg_name, "fontsize", 16);

title("log_{10} of the relative error for \\nu -> W(H^{(1)}_\\nu, H^{(2)}_\\nu)(\\nu w)", "fontsize", 16)

hold off;

%%%%%%%%%%%%%%%%%%%%
%%%% w -> J_nu(nu w)
figure(3);
hold on;

[_ i0] = min(abs(nu .- 1));
semilogy(r, err(:, i0), "+--");

[_ i1] = min(abs(nu .- 10));
semilogy(r, err(:, i1), "+--");

[_ i2] = min(abs(nu .- 100));
semilogy(r, err(:, i2), "+--");

xlabel("r", "fontsize", 16)
ylabel("relative error", "fontsize", 16)
axis tight;
grid on;
legend({[
  "\\nu = ", num2str(nu(i0))],
  ["\\nu = ", num2str(nu(i1))],
  ["\\nu = ", num2str(nu(i2))]
}, "fontsize", 16)

title("log_{10} of the relative error for w -> W(H^{(1)}_\\nu, H^{(2)}_\\nu)(\\nu w)", "fontsize", 16)

hold off;
