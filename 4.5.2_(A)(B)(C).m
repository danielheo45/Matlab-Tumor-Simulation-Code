clc; clear;

N1_0 = 0;
N2_0 = 5;
T_0 = 0.05;
G_0 = 0;

t_final = 24*30;
tspan = 0:0.01:t_final;
Dt = tspan;
D = zeros(1, length(Dt));

tic;
for i = 1:length(Dt)
    for j = 0:3
        for k = 0:4
            if Dt(i) >= 7*24*j+24*k && Dt(i) <= 7*24*j+24*k+24*0.3
                D(i) = 1;
            end
        end
    end
end

[t, y] = ode45(@TANs, tspan, [N1_0, N2_0, T_0, G_0],[],Dt,D);
N1_final = y(:, 1);
N2_final = y(:, 2);
T_final = y(:, 3);
G_final = y(:, 4);
toc;
%%
figure(1)
set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2)
plot(t./24, D, 'r', 'linewidth', 3)
xlabel('Time (days)')
ylabel('Radiation Dose (D)')

figure(2)
set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2)
plot(t./24, N1_final, 'b', t./24, N2_final, 'r', 'linewidth', 3)
xlabel('Time (days)')
ylabel('N2,N1 TANs')
legend({'N1 TANs', 'N2 TANs'})

figure(3)
set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2)
plot(t./24, T_final, 'r', 'linewidth', 3)
xlabel('Time (days)')
ylabel('Tumor population')

function dy = TANs(t, y, Dt, D)

mu_1 = 0.0289;
mu_2 = 0.0289;
mu_T = 0.005;
mu_G = 0.0289;

lambda_S = 0.0289;
lambda_6 = 0.000289;
lambda_G = 0.0289;
lambda_T = 0.0289;

k_1 = 0.1156;
k_2 = 1;
k_3 = 0.1156;
k_4 = 1;

alpha = 1;
beta = 1.5;
gamma = 0.015;
gamma_L = 2;

S = 0.2;
L = 0;

T_0 = 1;

alpha_RT = 0.03;
beta_RT = 0.003;

epsilon_N1 = 0.3;
epsilon_N2 = 0.3;

lambda_RT = 0.75;

D = interp1(Dt, D, t);

mu_RT = 1-exp(-alpha_RT*D-beta_RT*D^2);

N1 = y(1);
N2 = y(2);
T = y(3);
G = y(4);

dy = [lambda_S*S+k_1*k_2^2/(k_2^2+alpha*N2^2)-mu_1*N1-epsilon_N1*mu_RT*N1;
    lambda_6+lambda_G*G+k_3*k_4^2/(k_4^2+beta*N1^2)-mu_2*N2-epsilon_N2*mu_RT*N2;
    gamma*T*(1-T/T_0)-mu_T*N1*T-mu_RT*T;
    lambda_T*T-mu_G*G-gamma_L*L*G];

end