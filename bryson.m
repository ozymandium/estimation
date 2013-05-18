function [Ad,Qd] = bryson(Ac,Qc,Bw,Ts)
% Bryons's Trick - Convert the A and Q matrices from continuous to discrete
% 
% [Ad,Qd] = bryson(Ac,Qc,Bw)
% 
% INPUTS:
%       Ac = system dynamic matrix in continuous domain 
%       Bw = noise input matrix
%       Qc = process noise covariance matrix in continuous domain 
%       Ts = Time step (s)
% OUTPUTS:
%       Ad = system dynamic matrix in discrete domain 
%       Qd = process noise covariance matrix in discrete domain 

% number of states
nst = length(Ac);
if length(Bw)~=nst, error('must have noise input matrix length = # states'), end

S = [-Ac, Bw*Qc*Bw'; zeros(size(Ac)), Ac'];
C_bryson = expm(S.*Ts);

Ad = C_bryson( nst+1:2*nst, nst+1:2*nst )';
% % if Ad ~= 
% %     warning('Ad Bryson not equal to Ad');
% % end

Qd = Ad*C_bryson(1:nst,nst+1:2*nst);

end