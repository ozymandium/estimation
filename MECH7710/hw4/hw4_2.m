%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MECH 7710 HW 4, problem 2
%   using Michael Wooten's HW as a guide.
%   `aft' refers to occuring after the measurement update.
%   `bef' refers to ocurring before the measurement update.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
data=load('hw4_2.txt'); t=data(:,1); y=data(:,2); clear data

%% System/Setup
Ts = 0.1;
time = 0:Ts:150;
tlen = length(time);

Ac = 0; Ad = exp(Ac*Ts);
Bc = 0; Bd = exp(Bc*Ts);
Qd = 0.1; %! Tune this
Cc = 1; Cd = exp(Cc*Ts);
Rc = 1; Rd = exp(Rc*Ts);

v = 0 + 1.*randn(tlen,1); % generate Sensor noise
x = zeros(1,tlen); % state
y = zeros(1,tlen); % measurement

% system
for k = 1:tlen-1
    x(k+1) = Ad*x(k);
    y(k) = Cd*x(k) + Rd*v(k);
end


%% a

[L_ss,Pbef_ss,Paft_ss,poles] = dlqe(0,0,1,Qd,1);

% simulate
Pbef = zeros(1,tlen);
Pbef(1) = 1; % use steady state P value for initial P guess
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
L(k+1) = Paft(k+1)*Cd'/Rd;
x_aft(k+1) = x_bef(k) + L(k)*(y(k)-Cd*x_bef(k)); % update estimate

figa = namefig('a');
subplot(2,1,1); plot(time,L); grid on; title('Kalman Gain');
subplot(2,1,2); plot(time,x_aft); grid on; title('bias estimate');
