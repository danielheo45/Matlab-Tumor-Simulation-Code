clc; clear; close all;

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

G = 0.5;

h = 0.01;

for i = 1:501
    N1(i) = (i-1)*h;
    N2(i) = (lambda_6+lambda_G*G+k_3*k_4^2/(k_4^2+beta*N1(i)^2))/mu_2;
    S(i) = (mu_1*N1(i)-k_1*k_2^2/(k_2^2+alpha*N2(i)^2))/lambda_S;
end

hold on
figure(1)
plot(S, N1, S, N2, 'linewidth', 3)
xlabel('IFN-beta', 'fontsize', 30)
ylabel('Steady state', 'fontsize', 30)
legend({'N1 TANs', 'N2 TANs'}, 'fontsize', 30)
axis([0 0.6 0 5])
% plot(S, N2, 'b', 'linewidth', 3)
% xlabel('IFN-beta', 'fontsize', 30)
% ylabel('N2 Concentration', 'fontsize', 30)