clc; clear; close all;

% log(2)/24 = 0.0289
mu_1 = 0.0289;
mu_2 = 0.0289;

lambda_S = 0.0289;
lambda_6 = 0.000289;
lambda_G = 0.0289;

k_1 = 0.1156;
k_2 = 1;
k_3 = 0.1156;
k_4 = 1;

alpha = 1;
beta = 1.5;

S = 0;
G = 0.5;

% N2 nullcline
N1_base = 0:0.01:5;
N2 = (lambda_6+lambda_G*G+k_3*k_4^2./(k_4^2+beta*N1_base.^2))./mu_2;

% N1 nullcline
N2_base = 0:0.01:5;
N1 = (lambda_S*S+k_1*k_2^2./(k_2^2+alpha*N2_base.^2))./mu_1;

% [x, y] = meshgrid(linspace(0, 5, 100), linspace(0, 5, 100));
% 
% dx = lambda_S*S+k_1*k_2^2./(k_2^2+alpha*y.^2)-mu_1*x;
% dy = lambda_6+lambda_G*G+k_3*k_4^2./(k_4^2+beta*x.^2)-mu_2*y;
% 
% norm = sqrt(dx.^2+dy.^2);
% dx_norm = dx./norm;
% dy_norm = dy./norm;

figure(1)
% quiver(x, y, dx_norm, dy_norm)
% hold on
plot(N1_base, N2, 'r', N1, N2_base, 'b', 'linewidth', 3)
xlabel('N1 TANs', 'fontsize', 30)
ylabel('N2 TANs', 'fontsize', 30)
legend({'N2 nullcline', 'N1 nullcline'}, 'fontsize', 30)
