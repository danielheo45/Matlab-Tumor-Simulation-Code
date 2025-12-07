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

uS_total = 1.6;
uL = 0;

period = 2:1:15;

tic;
% every 5 day
% for i=1:length(St)
%     if (St(i)>=0 && St(i)<1*24)
%         injS(i) = uS;
%     elseif (St(i)>=5*24 && St(i)<6*24)
%         injS(i) = uS;
%     elseif (St(i)>=10*24 && St(i)<11*24)
%         injS(i) = uS;
%     elseif (St(i)>=15*24 && St(i)<16*24)
%         injS(i) = uS;
%     elseif (St(i)>=20*24 && St(i)<21*24)
%         injS(i) = uS;
%     elseif (St(i)>=25*24 && St(i)<26*24)
%         injS(i) = uS;
%     else
%         injS(i) = 0;
%     end
% end
% for i=1:length(Lt)
%     injL(i) = 0;
% end
% [t, y] = ode45(@TANs, tspan, [N1_0, N2_0, T_0, G_0, L_0, S_0],[],St,injS,Lt,injL);

% every 3 day
% for i=1:length(St)
%     if (St(i)>=0 && St(i)<1*24)
%         injS(i) = uS;
%     elseif (St(i)>=3*24 && St(i)<4*24)
%         injS(i) = uS;
%     elseif (St(i)>=6*24 && St(i)<7*24)
%         injS(i) = uS;
%     elseif (St(i)>=9*24 && St(i)<10*24)
%         injS(i) = uS;
%     elseif (St(i)>=12*24 && St(i)<13*24)
%         injS(i) = uS;
%     elseif (St(i)>=15*24 && St(i)<16*24)
%         injS(i) = uS;
%     elseif (St(i)>=18*24 && St(i)<19*24)
%         injS(i) = uS;
%     elseif (St(i)>=21*24 && St(i)<22*24)
%         injS(i) = uS;
%     elseif (St(i)>=24*24 && St(i)<25*24)
%         injS(i) = uS;
%     elseif (St(i)>=27*24 && St(i)<28*24)
%         injS(i) = uS;
%     else
%         injS(i) = 0;
%     end
% end
for j = 1:length(period)
    for k = 1:length(St)
        if mod(St(k), 24*period(j)) >= 0 && mod(St(k), 24*period(j)) <= 24
            injS(k) = uS_total/ceil(30/period(j));
        else
            injS(k) = 0;
        end
    end
    for k = 1:length(Lt)
        injL(k) = 0;
    end
    [t, y] = ode45(@TANs, tspan, [N1_0, N2_0, T_0, G_0, L_0, S_0],[],St,injS,Lt,injL);
    T_data(j) = y(end, 3);
end
toc;
%%
figure(1)
hold on
plot3(uS_total./ceil(30./period), period, 3*ones(length(period)), 'linewidth', 3)
% 
% figure(2)
% set(gcf,'Color',[1,1,1]);
% set(0,'DefaultAxesFontName', 'Times New Roman')
% set(0,'DefaultAxesFontSize', 25)
% set(0,'DefaultTextFontname', 'Times New Roman')
% set(0,'DefaultTextFontSize', 25)
% set(0,'defaultaxeslinewidth',2)
% plot(t./24,y(:,1),'b',t./24,y(:,2),'r','linewidth',3)
% xlabel('Time (Days)','fontsize',35)
% ylabel('TANs','fontsize',35)
% legend('N1 TANs','N2 TANs')
% 
% figure(3)
% set(gcf,'Color',[1,1,1]);
% set(0,'DefaultAxesFontName', 'Times New Roman')
% set(0,'DefaultAxesFontSize', 25)
% set(0,'DefaultTextFontname', 'Times New Roman')
% set(0,'DefaultTextFontSize', 25)
% set(0,'defaultaxeslinewidth',2)
% plot(t./24,y(:,3),'r','linewidth',3)
% xlabel('Time (Days)','fontsize',35)
% ylabel('Tumor','fontsize',35)

% [X, Y] = meshgrid(uS_total, period);

% figure(1)
% set(gcf,'Color',[1,1,1]);
% set(0,'DefaultAxesFontName', 'Times New Roman')
% set(0,'DefaultAxesFontSize', 25)
% set(0,'DefaultTextFontname', 'Times New Roman')
% set(0,'DefaultTextFontSize', 25)
% set(0,'defaultaxeslinewidth',2)
% surf(X, Y, T_data')
% axis([0.1 0.52 2 15])
% view(0, 90)
% colorbar
% xlabel('u_S')
% ylabel('period')

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