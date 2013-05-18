clear all; close all; clc
format compact
rng(1);

%% Constants
m = 10;
b = 7;
k = 10;

ts = 0.1;
time = 0:ts:20;
nt = length(time);
nsim = 1000;

%% System

Ac = [-b/m,-k/m;1,0];
Bc = [1/m;0];
Cc = [0,1];
Dc = 0;

[Ad,Bd,Cd,Dd] = c2dm(Ac,Bc,Cc,Dc, ts,'zoh');

plant_c = ss(Ac,Bc,Cc,Dc);
plant_d = c2d(plant_c, ts, 'zoh');
plant_tf = tf(1,[m,b,k]);
    
nst = length(Ac); % # states
nin = length(Dc); % # inputs
nms = min(size(Cd));

Bw = [0.5,0;0,0.5]; % noise input matrix
% Bwd= expm(Bwc.*ts);
Qc = cov(randn(2)); % process noise covariance in continuous

% % Bryson's Trick
S = [-Ac, Bw*Qc*Bw'; zeros(length(Ac)), Ac']; % Bryson's trick
C_bryson = expm(S.*ts);
Ad_bryson = C_bryson( nst+1:2*nst, nst+1:2*nst )';
if Ad_bryson ~= Ad, warning('Ad Bryson not equal to Ad'); end
Qd = Ad_bryson*C_bryson(1:nst,nst+1:2*nst); % process noise covariance in discrete; assume constant

% Qd_actual = zeros(nsim,1);
% for n = 1:nsim
%     Qd_actual(n) = 
% end

Rc = 10; % measurement noise covariance in continuous
Rd = exp(Rc*ts); % measurement noise covariance in discrete

% % estimate Pss
% Pss = lyap(Ac,Qc);
% Pssd = dlyap(Ad,Qd);

%% Simulation

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
% % should be able to start P from anywhere and have it converge
for n = 1:nsim
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

% % Mean over all runs
x_mean = mean(x,3);
y_mean = mean(y,3);
x_aft_mean = mean(x_aft,3);
x_err_mean = mean(x_err,3);

% % Ideal
[y_ideal, t_ideal, x_ideal] = lsim(plant_d, u(:,:,1), time, x(:,1,1));


%% Plotting

fig_sim = figure('Name','Simulation');
subplot(2,1,1)
plot(time,x_mean(2,1:end-1),...
     time,y_mean(1,1:end-1),...
     time,x_aft_mean(2,:))
legend('State','Meas','Est'); grid on; title('Filter')
subplot(2,1,2)
plot(time,x_err_mean)
title('Error Estimate'); grid on

fig_ideal = figure('Name','Ideal');
plot(t_ideal,y_ideal,...
     t_ideal,x_ideal(:,2))
legend('Meas','State'); grid on; title('Ideal')







