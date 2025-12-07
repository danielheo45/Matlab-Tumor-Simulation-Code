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

epsilon_N1 = 0.3;
epsilon_N2 = 0.3;

alpha_RT = 0.03;
beta_RT = 0.003;
D = [0 1 3 5];

for j = 1:4
    mu_RT = 1-exp(-alpha_RT*D(j)-beta_RT*D(j)^2);
    
    for i = 1:501
        N1(i) = (i-1)*h;
        N2(i) = (lambda_6+lambda_G*G+k_3*k_4^2/(k_4^2+beta*N1(i)^2))/(mu_2+epsilon_N2*mu_RT);
        S(i) = (mu_1*N1(i)-k_1*k_2^2/(k_2^2+alpha*N2(i)^2)+epsilon_N1*mu_RT*N1(i))/lambda_S;
    end
    N1_final(:, j) = N1;
    N2_final(:, j) = N2;
    S_final(:, j) = S;
end

hold on
D_x = 0:0.01:5;
mu_RT_y = 1.-exp(-alpha_RT.*D_x-beta_RT.*D_x.^2);
figure(1)
plot([0 5],[mu_RT_y(1) mu_RT_y(1)], [0 5],[mu_RT_y(101) mu_RT_y(101)], [0 5],[mu_RT_y(301) mu_RT_y(301)], [0 5],[mu_RT_y(501) mu_RT_y(501)], 'linewidth', 3)
hold on
plot(D_x, mu_RT_y, 'r', 'linewidth', 3)
xlabel('Radiation Dose (D)', 'fontsize', 25)
ylabel('$\mu_{RT}$', 'Interpreter', 'latex', 'fontsize', 30)
axis([0 5 -0.01 0.25])

figure(2)
plot(S_final(:, 1), N1_final(:, 1), S_final(:, 2), N1_final(:, 2), S_final(:, 3), N1_final(:, 3), S_final(:, 4), N1_final(:, 4), 'linewidth', 3)
xlabel('IFN-beta', 'fontsize', 25)
ylabel('Steady state (N1)', 'fontsize', 25)
legend({'D=0', 'D=1', 'D=3', 'D=5'}, 'fontsize', 25)
axis([0 0.6 0 5])

figure(3)
plot(S_final(:, 1), N2_final(:, 1), S_final(:, 2), N2_final(:, 2), S_final(:, 3), N2_final(:, 3), S_final(:, 4), N2_final(:, 4), 'linewidth', 3)
xlabel('IFN-beta', 'fontsize', 25)
ylabel('Steady state (N2)', 'fontsize', 25)
legend({'D=0', 'D=1', 'D=3', 'D=5'}, 'fontsize', 25)
axis([0 0.6 0 5])