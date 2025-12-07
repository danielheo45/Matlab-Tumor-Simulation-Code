clc; clear;

N1_0 = 0;
N2_0 = 5;
T_0 = 0.05;
G_0 = 0;
L_0 = 0;
S_0 = 0;

t_final = 24*30;
tspan = 0:0.01:t_final;
St = tspan;
Lt = tspan;

uS = 0.3;
uL = 0.1;

p = unique(perms([0 0 0 1 1 1]), 'rows');
p_label = unique(perms(['L' 'L' 'L' 'S' 'S' 'S']), 'rows');

tic;
for i = 1:20
    for j = 1:length(St)
        if St(j) >= 24*0 && St(j) <= 24*1
            injS(i, j) = uS*p(i, 1);
            injL(i, j) = uL*(1-p(i, 1));
        elseif St(j) >= 24*5 && St(j) <= 24*6
            injS(i, j) = uS*p(i, 2);
            injL(i, j) = uL*(1-p(i, 2));
        elseif St(j) >= 24*10 && St(j) <= 24*11
            injS(i, j) = uS*p(i, 3);
            injL(i, j) = uL*(1-p(i, 3));
        elseif St(j) >= 24*15 && St(j) <= 24*16
            injS(i, j) = uS*p(i, 4);
            injL(i, j) = uL*(1-p(i, 4));
        elseif St(j) >= 24*20 && St(j) <= 24*21
            injS(i, j) = uS*p(i, 5);
            injL(i, j) = uL*(1-p(i, 5));
        elseif St(j) >= 24*25 && St(j) <= 24*26
            injS(i, j) = uS*p(i, 6);
            injL(i, j) = uL*(1-p(i, 6));
        else
            injS(i, j) = 0;
            injL(i, j) = 0;
        end
    end
    [t, y] = ode45(@TANs, tspan, [N1_0, N2_0, T_0, G_0, L_0, S_0],[],St,injS(i, :),Lt,injL(i, :));
    N1_final(i, :) = y(:, 1);
    N2_final(i, :) = y(:, 2);
    T_final(i, :) = y(:, 3);
    G_final(i, :) = y(:, 4);
    L_final(i, :) = y(:, 5);
    S_final(i, :) = y(:, 6);
end
toc;
%%
% for i=1:20
%     subplot(4,5,i)
%     plot(t, S_final(i, :), 'r', t, G_final(i, :), 'b', 'linewidth', 3)
% end

figure(1)
set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2)
plot(t./24, injS(16, :), 'b', t./24, injL(16, :), 'r', 'linewidth', 3)
xlabel('Time(days)')
ylabel('injection rate(u_S,u_L)')
legend({'u_S', 'u_L'})
title('SLSSLL')

figure(2)
set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2)
plot(t./24, injS(20, :), 'b', t./24, injL(20, :), 'r', 'linewidth', 3)
xlabel('Time(days)')
ylabel('injection rate(u_S,u_L)')
legend({'u_S', 'u_L'})
title('SSSLLL')

figure(3)
set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2)
plot(t./24, N1_final(16, :), 'b', t./24, N2_final(16, :), 'r', 'linewidth', 3)
xlabel('Time(days)')
ylabel('N2,N1 TANs')
legend({'N1 TANs', 'N2 TANs'})
title('SLSSLL')

figure(4)
set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2)
plot(t./24, N1_final(20, :), 'b', t./24, N2_final(20, :), 'r', 'linewidth', 3)
xlabel('Time(days)')
ylabel('N2,N1 TANs')
legend({'N1 TANs', 'N2 TANs'})
title('SSSLLL')

figure(5)
set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2)
plot(t./24, T_final(16, :), 'b', t./24, T_final(20, :), 'r', 'linewidth', 3)
xlabel('Time(days)')
ylabel('Tumor population')
legend({'SLSSLL', 'SSSLLL'})

figure(6)
rmax = max(T_final(:, end));
theta = linspace(0, 2*pi);
polarplot(theta, rmax*ones(size(theta)), 'k', 'linewidth', 4)
hold on
polarplot([T_final(:, end)' T_final(1, end)], 'r', 'linewidth', 2)
thetaticks(0:18:342)
rlim([0 rmax])
rticklabels([])
rticks(0:rmax/5:rmax)
thetaticklabels(p_label)
set(gca, 'FontSize', 15, 'FontName', 'Times New Roman')
title('Tumor population', 'FontSize', 25)

function dy = TANs(t, y,St,injS,Lt,injL)

mu_1 = 0.0289;
mu_2 = 0.0289;
mu_T = 0.005;
mu_G = 0.0289;
mu_L = 0.231;
mu_S = 0.1386;

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

T_0 = 1;

injS = interp1(St,injS,t);
injL = interp1(Lt,injL,t);

N1 = y(1);
N2 = y(2);
T = y(3);
G = y(4);
L = y(5);
S = y(6);

dy = [lambda_S*S+k_1*k_2^2/(k_2^2+alpha*N2^2)-mu_1*N1;
    lambda_6+lambda_G*G+k_3*k_4^2/(k_4^2+beta*N1^2)-mu_2*N2;
    gamma*T*(1-T/T_0)-mu_T*N1*T;
    lambda_T*T-mu_G*G-gamma_L*L*G;
    injL-mu_L*L;
    injS-mu_S*S];

end