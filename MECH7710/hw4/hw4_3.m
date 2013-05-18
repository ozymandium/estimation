%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MECH 7710 HW 4, problem 3
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

% read data
data_raw = load('./hw4_3.txt');
time=data_raw(:,1); east=data_raw(:,2); north=data_raw(:,3); psi=data_raw(:,4); gyro=data_raw(:,5); radar=data_raw(:,6); clear data_raw

tlen = length(time);
Ts = 0.2;
nst = 5;

%% Setup

x = zeros(nst,tlen);
y = zeros(nst,tlen);
Ac = 

%% Run
for k=1:tlen
    
    
    
end