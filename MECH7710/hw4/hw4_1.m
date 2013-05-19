%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MECH 7710 HW 4, problem 1
%   using Michael Wooten's HW as a guide.
%   `aft' refers to occuring after the measurement update.
%   `bef' refers to ocurring before the measurement update.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%% System/Setup
Ts = 0.1;

Ac = [0,1;-1,-1.4];
Gc = [0;1]; % process noise input matrix
Gd = expm(Gc.*Ts); 
Qc = 5^2;
[Ad,Qd] = bryson(Ac,Qc,Gc,Ts); Qd

Cc = [1 0];
Cd = exp(Cc*Ts);
Rc = 1^2;
Rd = expm(Rc*Ts)

%% a - Simulate Controlled system
hw4
time = 0:Ts:100;
tlen = length(time);

v = 0 + 1.*randn(tlen,1); % generate Sensor noise
w = 0 + 2.*randn(tlen,1); % generate Process noise
x = zeros(2,tlen); % state
y = zeros(1,tlen); % measurement

for k = 1:tlen-1
    x(:,k+1) = Ad*x(:,k) + Gd*w(k);
    y(:,k) = Cd*x(:,k) + v(k);
end

%% b - find true Q,Qd,Rd

Rd_true = cov(v);
Rc_true = log(Rd_true)/Ts;
fprintf('True Rd:  %f\n',Rd_true)
Qd_true = cov(w);
Qc_true = log(Qd_true)/Ts;
fprintf('True '),Qc_true, fprintf('True '), Qd_true

%% c - Steady State calculation

[L_ss,Pbef_ss,Paft_ss,poles] = dlqe(Ad,eye(2),Cd,Qd,Rd)

%% d - Steady Stahw4te Kalman Filter

Pbef = zeros(2,2,tlen);
Pbef(:,:,1) = Pbef_ss; % use steady state P value for initial P guess
Paft = zeros(2,2,tlen);
L = zeros(2,tlen);
x_bef = zeros(2,tlen);
x_aft = zeros(2,tlen);

% Run estimation for comparison
for k=1:tlen-1
    Paft(:,:,k) = inv( inv(Pbef(:,:,k)) + Cd'/Rd*Cd ); % Estimate Uncertainty
    % usually L(:,k) = Paft(:,:,k)*Cd'/Rd; % Kalman Gain
    L(:,k) = L_ss;
    x_aft(:,k) = x_bef(:,k) + L(:,k)*(y(:,k)-Cd*x_bef(:,k)); % update estimate
    
    % Propagate to next time
    x_bef(:,k+1) = Ad*x_aft(:,k);
    Pbef(:,:,k+1) = Ad*Paft(:,:,k)*Ad' + Qd;
end
x_aft(:,k+1) = x_bef(:,k) + L(:,k)*(y(:,k)-Cd*x_bef(:,k)); % update estimate

% norm of errors
N = sqrt( sum( std(x-x_aft,0,2).^2 ) )

figa = namefig('Comparison of state, estimate, meas');
subplot(3,1,1:2)
plot(time,x(1,:),...
     time,y,...
     time,x_aft(1,:),...
     'LineWidth',2)
legend('State','Meas','Est'); title('Comparison of state, estimate, meas'); xlabel('Time'); grid on; hold on    
subplot(3,1,3)
plot(time,x(1,:)-x_aft(1,:));
grid on; title('errors')

%% e - Qd/Rd ratio

% Larger Qd/Rd ratios give lower values of the the norm error. The
% difference stops being appreciable around Qc=5-10.












