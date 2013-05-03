%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Kalman Filter for Linear System - Mass,Spring,Damper
%   - comparison for ideal, unnoisy system with lsim
%   - should work for any linear system by changing constants 
%       Continuous state space matrices
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
format compact
rng(1);

%% Constants
m = 5;
b = 4;
k = 7;

ts = 0.05;
time = 0:ts:20;
nt = length(time);
nsim = 500;

%% System

Ac = [-b/m,-k/m;1,0];
Bc = [1/m;0];
Cc = [0,1];
Dc = 0;

[Ad,Bd,Cd,Dd] = c2dm(Ac,Bc,Cc,Dc, ts,'zoh');

plant_c = ss(Ac,Bc,Cc,Dc);
plant_d = c2d(plant_c, ts, 'zoh');
% plant_tf = tf(1,[m,b,k]);

nst = length(Ac); % # states
nin = length(Dc); % # inputs
nms = min(size(Cd));

Bw = [1,0;0,1]; % noise input matrix (continuous)
Qc = [1,0;0,1]; % process noise covariance in continuous

% % Bryson's Trick
S = [-Ac, Bw*Qc*Bw'; zeros(length(Ac)), Ac']; % Bryson's trick
C_bryson = expm(S.*ts);
Ad_bryson = C_bryson( nst+1:2*nst, nst+1:2*nst )';
if Ad_bryson ~= Ad, warning('Ad Bryson not equal to Ad'); end
Qd = Ad_bryson*C_bryson(1:nst,nst+1:2*nst); % process noise covariance in discrete; assume constant

Rc = 10; % measurement noise covariance in continuous
Rd = exp(Rc*ts); % measurement noise covariance in discrete

% % estimate Pss
% Pss = lyap(Ac,Qc);
% Pssd = dlyap(Ad,Qd);

%% Simulation - KF

% % Input
in_mag = 10;
u = in_mag*ones(nin,nt,nsim); % input - step

x = zeros(nst,nt,nsim); % states
x_bef = zeros(nst,nt,nsim); % state estimates before update
x_aft = zeros(nst,nt,nsim); % state estimates after update
x_err = zeros(nms,nt,nsim); % state estimate errors
y = zeros(nms,nt,nsim); % measurements
w = zeros(2,nt,nsim); % process noise
v = randn(nms,nt,nsim); % measurement noise
P_bef = zeros(nst,nst,nt,nsim);
P_aft = zeros(nst,nst,nt,nsim);
L = zeros(nst,nt,nsim);

% % Initial conditions
for n = 1:nsim % should be able to start P from anywhere and have it converge???
    P_bef(:,:,1,n) = dlyap(Ad,Qd); % start P at the steady state value
end

for n = 1:nsim
    for k = 1:nt       
        % time update
        w(:,k,n) = sqrtm(Qd)*randn(nst,1);
        x(:,k+1,n) = Ad*x(:,k,n) + Bd*u(1,k,n) + Bw*w(:,k,n);
        y(:,k+1,n) = Cd*x(:,k,n) + Rd*v(1,k,n);
        
        % Kalman gain
        P_aft(:,:,k,n) = inv( P_bef(:,:,k,n) + Cd'*inv(Rd)*Cd );
        L(:,k,n) = P_aft(:,:,k,n)*Cd'*inv(Rd);
        
        % measurement update
        x_err(:,k,n) = y(1,k,n) - Cd*x_bef(:,k,n);
        x_aft(:,k,n) = x_bef(:,k,n) + L(:,k,n)*x_err(:,k,n);
%         x_aft(:,k,n) = L(:,k,n)'*x_bef(:,k,n) + L(:,k,n)*y(:,k,n);

        % Next step prediction
        x_bef(:,k+1,n) = Ad*x_aft(:,k,n) + Bd*u(:,k,n);
        P_bef(:,:,k+1,n) = Ad*P_aft(:,:,k,n)*Ad' + Qd; % can even do this before running the simulation
    end        
end


%% Post Mortem

% % Mean over all runs
x_mean = mean(x,3);
y_mean = mean(y,3);
x_aft_mean = mean(x_aft,3);
x_err_mean = mean(x_err,3);
x_err_true = x_aft-x(:,1:end-1,:); % true error in estimate
x_err_true_mean = mean(x_err_true,3);

% % Ideal
[y_ideal, t_ideal, x_ideal] = lsim(plant_d, u(:,:,1), time, x(:,1,1));

%% Plotting

fig_sim = figure('Name','Simulation');
subplot(3,1,[1 2])
plot(time,x_mean(2,1:end-1),...
     time,y_mean(1,1:end-1),...
     time,x_aft_mean(2,:),...
     t_ideal,y_ideal,'k')
legend('State','Meas','Est','Ideal'); grid on; title('Filter')
subplot(3,1,3)
plot(time,x_err_mean,...
     time,x_err_true_mean(2,:))
title('Errors'); grid on; legend('Est.','True')




