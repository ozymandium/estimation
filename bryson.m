function [Ad,Qd] = bryson(Ac,Qc,G,Ts)
% Bryons's Trick - Convert the A and Q matrices from continuous to discrete
% 
% [Ad,Qd] = bryson(Ac,Qc,G,Ts)
% 
% INPUTS:
%       Ac = system dynamic matrix in continuous domain 
%       Qc = process noise covariance matrix in continuous domain 
%       G = noise input matrix
%       Ts = Time step (s)
% OUTPUTS:
%       Ad = system dynamic matrix in discrete domain 
%       Qd = process noise covariance matrix in discrete domain 

% number of states
n = length(Ac);
if length(G)~=n, error('must have noise input matrix length = # states'), end

S_bryson = [-Ac, G*Qc*G'; zeros(n), Ac'];
C_bryson = expm(S_bryson.*Ts);

Ad = C_bryson( n+1:2*n, n+1:2*n )';
if Ad~=exp(Ac*Ts)
    warning('bryson:Ad_Error','Ad Bryson not equal to exp(Ac*Ts)');
end

Qd = Ad*C_bryson(1:n,n+1:2*n);

end