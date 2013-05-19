%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Simple Kalman Filter for Mass-Spring-Damper System
%   `aft' refers to occuring after the measurement update (a posteriori).
%   `bef' refers to ocurring before the measurement update (a priori).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
% format compact
% rng(1);

%% Constants
mass = 10;
damping = 7;
spring = 5;

Ts = 0.05;
time = 0:Ts:20;
tlen = length(time);

in_mag = 20;
wstd = 0.05;
vstd = 0.1;

%% System

Ac = [-damping/mass, -spring/mass; 1,0];
Bc = [1/mass; 0];
G = eye(2,2); % noise input matrix
Cc = [0,1]; Cd = Cc;
Dc = 0; Dd = Dc;

n = length(Ac); % # states
l = length(Dc); % # inputs
m = size(Cc,1); % # measurements

% [Ad,Bd,Cd,Dd] = c2dm(Ac,Bc,Cc,Dc, Ts,'zoh');
plant_c = ss(Ac,Bc,Cc,Dc);
plant_d = c2d(plant_c, Ts, 'zoh');
    
% process noise covariance in continuous
% Qc = wstd^2*eye(n,n) % this causes estimate to follow ideal behavior of
% unnoisy system
Qc = 0.3*eye(n,n)

% measurement noise covariance in continuous-time
% Rc = vstd^2 % this causes estimate to follow ideal behavior of
% unnoisy system
Rc = 0.001
Rd = exp(Rc*Ts); % measurement noise covariance in discrete-time

[Ad,Qd] = bryson(Ac,Qc,G,Ts)
Bd = Ac\(Ad-eye(n,n))*Bc;


% Steady State
[Lss,Pbef_ss,Paft_ss,poles] = dlqe(Ad,eye(n,n),Cd,Qd,Rd);

%% Simulation

% % Inputs and Noises
u = in_mag*ones(l,tlen); % input - step
w = zeros(n,tlen); % process noise in discrete-time
v = zeros(m,tlen);  
w(2,:) = 0 + wstd.*randn(1,tlen); % only apply noise to position
v(1,:) = 0 + vstd.*randn(1,tlen); % measurement noise

% Preallocation
x = zeros(n,tlen); % states
y = zeros(m,tlen); % measurements
L = zeros(n,m,tlen);
S = zeros(m,tlen); % Residual/Innovation
xaft = zeros(n,tlen); % state estimates after update
Paft = zeros(n,n,tlen);
xbef = zeros(n,tlen); % state estimates before update
Pbef = zeros(n,n,tlen);

% Initial conditions
% % should be able to start P from anywhere and have it converge
Pbef(:,:,1) = Pbef_ss;

% plant simulation
for k=1:tlen
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k) + w(:,k);
    y(:,k+1) = Cd*x(:,k) + v(1,k);   
end

% Estimator Simulation
for k = 1:tlen         
    % Kalman gain
    L(:,:,k) = Pbef(:,:,k)*Cd'/(Cd*Pbef(:,:,k)*Cd'+Rd);

    % measurement update
    S(:,k) = y(:,k) - Cd*xbef(:,k);
    xaft(:,k) = xbef(:,k) + L(:,:,k)*S(:,k);
    Paft(:,:,k) = (eye(n,n)-L(:,:,k)*Cd)*Pbef(:,:,k);

    % Time update
    xbef(:,k+1) = Ad*xaft(:,k) + Bd*u(:,k);
    Pbef(:,:,k+1) = Ad*Paft(:,:,k)*Ad' + Qd; % can even do this before running the simulation
end        

%% Post Mortem

% Ideal
[y_ideal, t_ideal, x_ideal] = lsim(plant_d, u(:,:,1), time, x(:,1,1));

xaft_err = x(:,1:tlen) - xaft;
xbef_err = x(:,1:tlen) - xbef(:,1:tlen);
% Norm of std of errors
Naft = sqrt(sum(std(xaft_err,0,2).^2))
Nbef = sqrt(sum(std(xbef_err,0,2).^2))

%% Plotting

fig_sim = namefig('Simulation');
subplot(3,1,1:2)
plot(time,x(2,1:tlen),...
     time,y(1,1:tlen),...
     time,xbef(2,1:tlen),...
     time,xaft(2,1:tlen),...
     t_ideal,y_ideal,...
     'LineWidth',2)
legend('State','Meas','Est a pri','Est a post','Model'); grid on; title('Filter')
subplot(3,1,3)
plot(time,xbef_err(2,:),...
     time,xaft_err(2,:))
legend('a priori','a posteriori'); title('Error Estimate'); grid on

% fig_ideal = namefig('Ideal');
% plot(t_ideal,y_ideal,...
%      t_ideal,x_ideal(:,2))
% legend('Meas','State'); grid on; title('Ideal')
% 
% 





