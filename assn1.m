% Assi[gnment 1, Bevly's discretization homework
% Author: Robert Cofield
% Translated from python implementation 2012 - 10 - 12
%   Generic mass/spring/damper system
%   Stable and Underdamped
%   m*xDD + b*xD + k*x = Fin
clear all
close all
clc

% % Parameters
m = 2;
b = 2;
k = 4;
Fin = 1.0;
t_max = 10;
dt_ana = 0.1;
dt_eul = 0.05;
dt_trp = 0.05;
dt_dsc = 1;

sys = tf(1,[m,b,k]);
    sys.outputn = {'Position'};
    sys.outputu = {'m'};
    sys.inputn = {'Force'};
    sys.inputu = {'N'};
[numd, dend] = c2dm(sys.num{1}, sys.num{1}, dt_dsc);
sysd = tf(numd,dend,dt_dsc);
    sysd.outputn = {'Position'};
    sysd.outputu = {'m'};
    sysd.inputn = {'Force'};
    sysd.inputu = {'N'};
[A,B,C,D] = tf2ss(sys.num{1}, sys.den{1});
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D, dt_dsc, 'zoh');

% % General System stuff
% Auto-functions
fbode = namefig('Bode - Continuous & Discrete');
subplot(1,2,1), bode(sys),                 grid on
subplot(1,2,2), dbode(numd, dend, dt_dsc), grid on

floc = namefig('Poles & Zeros - Continuous & Discrete');
subplot(1,2,1), rlocus(sys), grid on
subplot(1,2,2), pzmap(sysd), grid on
grid on

x_ss = Fin/k;
zeta = b/(2*sqrt(m*k));
% Underdamped Stuff
wn = sqrt(k/m);
wd = wn*sqrt(1-zeta^2);
tau = 1/(zeta*wn);
% Continuous Eigenvalues
s1 = (-b - sqrt(b^2 - 4*m*k))/(2*m);
s2 = (-b + sqrt(b^2 - 4*m*k))/(2*m);
% Discrete Eigenvalues
z1 = exp(s1*dt_dsc);
z2 = exp(s2*dt_dsc);
% Real Eigenvalues
s = real(s1); z = real(z1);

% % Print info
fprintf('Step Info for System\n'), disp(stepinfo(sys))
fprintf('Steady State Resp.:    %f\n', x_ss)
fprintf('Damp. Ratio:           %f\n', zeta)
fprintf('Natural Frquency:      %f\n', wn)
fprintf('Damped Nat. Freq:      %f\n', wd)
fprintf('Time Constant:         %f\n', tau)
fprintf('\nEigenvalues')
fprintf('\n\tContinuous:\n\t\t%s\n\t\t%s\n', num2str(s1), num2str(s2))
fprintf('\n\tDiscrete:\n\t\t%s\n\t\t%s\n', num2str(z1), num2str(z2))


% % Simulations
% Analytical
t_ana = 0:dt_ana:t_max;
x_ana = zeros(1,length(t_ana));
C1 = -0.25;
for n = 1:length(t_ana)
    x_ana(n) = C1*cos(sqrt(7)/2*t_ana(n))*exp(-t_ana(n)/2) + 0.25;
end

% Continuous - Euler
t_eul = 0:dt_eul:t_max;
xDD_eul = zeros(1,length(t_eul));
xD_eul = zeros(1,length(t_eul));
x_eul = zeros(1,length(t_eul));
for n = 2:length(t_eul)
    xDD_eul(n) = (Fin - b*xD_eul(n-1) - k*x_eul(n-1))/m;
    xD_eul(n)  = xD_eul(n-1) + xDD_eul(n-1)*dt_eul;
    x_eul(n)   = x_eul(n-1)  + xD_eul(n-1) *dt_eul;
end

%  Numerical - Trapezoid
t_trp = 0:dt_trp:t_max;
xDD_trp = zeros(1,length(t_trp));
xD_trp = zeros(1,length(t_trp));
x_trp = zeros(1,length(t_trp));
for n = 2:length(t_eul)
    xDD_trp(n) = (Fin - b*xD_trp(n-1) - k*x_trp(n-1))/m;
    xD_trp(n) = xD_trp(n-1) + (xDD_trp(n)+xDD_trp(n-1))*dt_trp/2;
    x_trp(n)  = x_trp(n-1)  + (xD_trp(n)+ xD_trp(n-1))  *dt_trp/2;
end

% Discrete - State Space 
t_dsc = 0:dt_dsc:t_max;
X1 = zeros(1,length(t_dsc)+1);
X2 = zeros(1,length(t_dsc)+1);
u = Fin;
x_dsc_ss = zeros(1,length(t_dsc));
for k = 2:length(t_dsc)+1
    X_ = [X1(k-1); X2(k-1)];
    X = Ad*X_ + Bd*u;
    Y = Cd*X_ + Dd*u;
    X1(k) = X(1);
    X2(k) = X(2);
    x_dsc_ss(k-1) = Y(1);
end

% Discrete - Euler
x1_dsc_eul = zeros(1, length(t_dsc));
x2_dsc_eul = x1_dsc_eul;
for k = 3:length(x1_dsc_eul)+2
    x1_dsc_eul(k) = u/2 - x1_dsc_eul(k-1) - 2*x1_dsc_eul(k-2);
    x2_dsc_eul(k) = u/2 - x2_dsc_eul(k-1) - 2*x2_dsc_eul(k-2);

    x1_dsc_eul(k-1) = x1_dsc_eul(k)*z1;
    x2_dsc_eul(k-1) = x2_dsc_eul(k)*z2;

    x1_dsc_eul(k-2) = x1_dsc_eul(k-1)*z1;
    x2_dsc_eul(k-2) = x2_dsc_eul(k-1)*z2;
end
figure
plot(t_dsc, x1_dsc_eul(1:end-2), t_dsc, x2_dsc_eul(1:end-2))

% % Auto-Simulations
fig_auto = namefig('AutoSims');
hold on
step(sys)
dstep(numd, dend, length(t_dsc))
%     t_dstp = max(t_dsc)*sort(st_dstp(:,1))';

% % Plotting
close(floc,fbode)
fig_comp = namefig('Step Response Simulations');
ax = axes('Color', [0.5, 0.5, 0.5]);
hold on
plot(t_ana, x_ana,    'm',...
     t_eul, x_eul,    'c',...
     t_trp, x_trp,    'k')
stairs(t_dsc, x_dsc_ss,'r')
% stairs(t_dstp, x_dstp, 'b')
grid on
title('Step Response')
xlabel('Time (s)')
ylabel('Position (m)')
legend(['Analytical: dt=',num2str(dt_ana)],...
       ['Cont Euler: dt=',num2str(dt_eul)],...
       ['Cont Trapz: dt=',num2str(dt_trp)],...
       ['Disc St Sp: dt=',num2str(dt_dsc)],...
       'Location','SouthEast')

% Discrete - Euler - Noise response






