%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MECH 7710 HW 4, problem 2
%   using Michael Wooten's HW as a guide.
%   `aft' refers to occuring after the measurement update.
%   `bef' refers to ocurring before the measurement update.
%   We know the true state is zero throughout the data set
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
data=load('hw4_2.txt'); time=data(:,1); y=data(:,2); clear data

%% System/Setup
Ts = 0.1;
tlen = length(time);

Ac = 0; Ad = exp(Ac*Ts);
Gc = 0; Gd = exp(Gc*Ts);
Qd = 0.01 %! Tune this

Cc = 1; Cd = exp(Cc*Ts);
Rc = 2; Rd = exp(Rc*Ts)

%% a/b

[L_ss,Pbef_ss,Paft_ss,poles] = dlqe(Ad,Gd,Cd,Qd,Rd)

% simulate
Pbef = zeros(1,tlen);
Pbef(1) = 10; %! tune this as well???
Paft = zeros(1,tlen);
L = zeros(1,tlen);
x_bef = ones(1,tlen);
x_aft = ones(1,tlen);

% Run estimation for comparison
for k=1:tlen-1
    Paft(k) = inv( inv(Pbef(k)) + Cd'/Rd*Cd ); % Estimate Uncertainty
    L(k) = Paft(k)*Cd'/Rd; % Kalman Gain
    % L(k) = L_ss;
    x_aft(k) = x_bef(k) + L(k)*(y(k)-Cd*x_bef(k)); % update estimate
    
    % Propagate to next time
    x_bef(k+1) = Ad*x_aft(k);
    Pbef(k+1) = Ad*Paft(k)*Ad' + Qd;
end
% last time step
Paft(k+1) = inv( inv(Pbef(k+1)) + Cd'/Rd*Cd ); % Estimate Uncertainty
L(k+1) = Paft(k+1)*Cd'/Rd;
x_aft(k+1) = x_bef(k) + L(k)*(y(k)-Cd*x_bef(k)); % update estimate

%% c/d

numd_filter = sqrt(Qd);
dend_filter = [1, -(1-sqrt(Qd))];
y0 = y(1)
yf = filter(numd_filter, dend_filter, y, y0);

%% plot

figa = namefig('a/b/c');
subplot(2,1,1); plot(time,L); grid on; title('Kalman Gain');
subplot(2,1,2)
plot(time,y,...
    time,x_aft,...
    time,yf,...
    'LineWidth',2)
grid on; title('bias estimate & meas vs filt meas'); legend('meas','est','filt'); hold on







