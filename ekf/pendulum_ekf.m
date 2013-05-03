%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Extended Kalman Filter for NonLinear System 
%   - Pendulum with Damper, input Torque
%   - Qc is estimated by cov(randn(< # states >))
%   - Comparison for actual nonlinear via trapezoid integration
%       - changing system means needing to change the nonlinear simulation    
%   - !!! need to propagate estimates of A,C,R !!!
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
format compact
% rng(1);
%% Parameters

m = 2;
r = 1;
b = 4;
g = 9.81;
J = m*r^2;

ts = 0.01;
reinit_P_int = 5; % how many time steps to go before reinitializing P
time = 0:ts:5;
nt = length(time);
nsim = 100;

in_mag = 10; % magnitude of input

%% System

Ac = [-b/J,-m*g/J;1,0];
Bc = [1/J;0];
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
S = [-Ac, Bw*Qc*Bw'; zeros(nst), Ac']; % Bryson's trick
C_bryson = expm(S.*ts);
Ad_bryson = C_bryson( nst+1:2*nst, nst+1:2*nst )';
if Ad_bryson ~= Ad, warning('Ad Bryson not equal to Ad'); end
Qd = Ad_bryson*C_bryson(1:nst,nst+1:end); % process noise covariance in discrete; assume constant

Rc = 1; % measurement noise covariance in continuous
% Rd = exp(Rc*ts); % measurement noise covariance in discrete
Rd = 1;
%% Simulation

% % Input
u = in_mag*ones(nin,nt,nsim); % input - step

x = zeros(nst,nt,nsim); % states
x_bef = zeros(nst,nt,nsim); % state estimates before update
x_aft = zeros(nst,nt,nsim); % state estimates after update
x_err = zeros(nms,nt,nsim); % state estimate errors
y = zeros(nms,nt,nsim); % measurements
w = zeros(nst,nt,nsim); % process noise
v = randn(nms,nt,nsim); % measurement noise
P_bef = zeros(nst,nst,nt,nsim);
P_aft = zeros(nst,nst,nt,nsim);
L = zeros(nst,nt,nsim);

% % Initial conditions
for n = 1:nsim % !? should be able to start P from anywhere and have it converge, but can't here...
    P_bef(:,:,1,n) = dlyap(Ad,Qd); % start P at the steady state value
    P_aft(:,:,1,n) = dlyap(Ad,Qd);
end

% % Make matrix A 3D so it can be updated in time
Ad_ = Ad;
Ad = zeros(nst,nst,nt,nsim);
for n=1:nsim
    Ad(:,:,1,n) = Ad_;
end

% % Actual simulation part
for n = 1:nsim
    for k = 2:nt       
        % Update state matrices - system-specific
        Ad(:,:,k,n) = Ad(:,:,k-1,n);
%         if x_aft(2,k,n) % update Ad index with sine if non-zero state
%             Ad(1,2,k,n) = Ad(1,2,k,n)*sin(x_aft(2,k,n))/x_aft(2,k,n); % propagate sine(theta) forward
%         end
        
        % time update for states
        w(:,k-1,n) = sqrtm(Qd)*randn(nst,1);
        x(:,k,n) = Ad(:,:,k)*x(:,k-1,n) + Bd*u(1,k-1,n);% + w(:,k-1,n);
        y(:,k,n) = Cd*x(:,k-1,n);% + Rd*v(1,k-1,n);
              
        % Time update
        P_bef_dot = Ad(:,:,k,n)*P_bef(:,:,k-1,n) + P_bef(:,:,k-1,n)*Ad(:,:,k,n)' + Qd; % !? should PCR-1CP be here?
        P_bef(:,:,k,n) = P_aft(:,:,k-1,n) + P_bef_dot*ts;
        x_bef_dot = Ad(:,:,k,n)*x_bef(:,k-1,n) + Bd*u(:,k-1,n);
        x_bef(:,k,n) = x_aft(:,k-1,n) + x_bef_dot*ts;

        % Kalman Gain
        L(:,k,n) = P_bef(:,:,k,n)*Cd'*inv(Cd*P_bef(:,:,k,n)*Cd'+Rd);

        % Measurement update
        x_err(:,k,n) = y(:,k,n) - Cd*x_bef(:,k,n); % estimate of state error 
        x_aft(:,k,n) = x_bef(:,k,n) + L(:,k,n)*x_err(:,k,n);
        if floor(k/reinit_P_int)==k 
            P_aft(:,:,k,n) = 0.5*(P_aft(:,:,k,n)+P_aft(:,:,k,n)');
        else
            P_aft(:,:,k,n) = ( eye(nst) - L(:,k,n)*Cd )*P_bef(:,:,k,n); 
        end
    end        
end

%% Simulation: Nonlinear system

thetaDD_nl = zeros(nsim,nt);
thetaD_nl = zeros(nsim,nt);
theta_nl = zeros(nsim,nt);

for n = 1:nsim
    for k = 1:nt-1 % use trapezoid integration
        thetaDD_nl(n,k+1) = -b/J*thetaD_nl(n,k) - m*g/J*sin(theta_nl(n,k)) + 1/J*u(1,k,n) + w(1,k,n);
        thetaD_nl(n,k+1) = thetaD_nl(n,k) + 0.5*sum(thetaDD_nl(n,k:k+1))*ts + w(2,k,n);
        theta_nl(n,k+1) = theta_nl(n,k) + 0.5*sum(thetaD_nl(n,k:k+1))*ts;
    end
end

%% Post Mortem

% % Mean over all runs
x_mean = mean(x,3);
y_mean = mean(y,3);
x_bef_mean = mean(x_bef,3);
x_aft_mean = mean(x_aft,3);
x_err_mean = mean(x_err,3);
x_err_true = x-x_aft; % true error in estimate
x_err_true_mean = mean(x_err_true,3);
theta_nl_mean = mean(theta_nl,1);

% % Ideal
[y_ideal, t_ideal, x_ideal] = lsim(plant_d, u(:,:,1), time, x(:,1,1));

%% Plotting

fig_sim = figure('Name','Simulation');
subplot(3,1,[1 2])
plot(time,x_mean(2,:),.....
     time,y_mean(1,:),...
     time,x_bef_mean(2,:),...
     t_ideal,y_ideal,'k',...
     time,theta_nl_mean,...
     'LineWidth',2)
legend('State','Meas','Est','Ideal Lin','Nonlinear', 'Location','Best'); grid on; title('Filter')
subplot(3,1,3)
plot(time,x_err_mean,...
     time,x_err_true_mean(2,:))
title('Errors'); grid on; legend('Est.','True')




