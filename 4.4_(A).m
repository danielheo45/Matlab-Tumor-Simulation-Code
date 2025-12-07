clc; clear; close all;

N1_0 = 0;
N2_0 = 5;
S_0 = 0;

t_final = 100000;

[t y] = ode45(@TANs, [0:0.01:t_final], [N1_0, N2_0, S_0]);

figure(1)
plot(t, y(:, 1), 'b', t, y(:, 2), 'r', t, y(:, 3), '--w', 'linewidth', 3)
xlabel('Time', 'fontsize', 30)
ylabel('N2, N1 steady states', 'fontsize', 30)
legend({'N1 TANs', 'N2 TANs', 'IFNbeta'})

function dy = TANs(t, y)

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
% S = 0.2;

N1 = y(1);
N2 = y(2);
S = y(3);

dy = [lambda_S*S+k_1*k_2^2./(k_2^2+alpha*N2.^2)-mu_1*N1;
    lambda_6+lambda_G*G+k_3*k_4^2./(k_4^2+beta*N1.^2)-mu_2*N2;
    1/100000];

end