% --------------------------------------------------------
% The follow codes computes the algebraic solution for X and T
% It also compute the interpolating values for them 
% Last change: August 15, 2019
% --------------------------------------------------------
clear
tic
% ----------------------------
% Parameters
% ----------------------------
M= 8; 
p = 10327 ; q = 27993 ;% for p/q = 1/3 + 1/31 + 1/301  
p = 18209 ; q = 65764 ; 
p = 6004; q = 18003; % p/q = 1/3 + 1/3001
p = 8001 ; q = 31988 ; % p / q = 1/4 + 1/7997
p = 18209 ; q = 65764 ; % p / q = 1/4 + 1/41 + 1/401

MM = 8 ;
L = 4.8 ; 
l = L / MM ; 
% l = 0.2 ; L = l*M ;

disc = 2^1; 
% Warning: Do not decrease N values otherwise X(0,t) will not be captured
% correctly 
N = 2^8*M ; 
% ----------------------------
T_run = l^2 / (2*pi) ; dt = T_run / (q+1);
c0 = sqrt(2*log(cosh(l/2))/pi) ;
cM=2*pi*c0^2 / (l * sqrt(1-exp(-pi*c0^2))) ;
s =linspace(-L/2,L/2,N+1);
t = linspace(0,T_run,q+1); 
rhoq = acosh(2*cosh(l/2)^(2/q)-1) ;
% X1fu_alg = zeros(N+1,q+1); X2fu_alg = zeros(N+1,q+1); X3fu_alg = zeros(N+1,q+1); 
% T1fu_alg = zeros(N+1,q+1); T2fu_alg = zeros(N+1,q+1); T3fu_alg = zeros(N+1,q+1); 
X0 = zeros(q+1,3) ; mtmXX = zeros(q+1,3); mtmXXq = zeros(q+1,3); 
mtmTT = zeros(q+1,3); mtmTTq = zeros(q+1,3); cpT = zeros(q+1,3); mtmTTM  = zeros(q+1,3);

% ----------------------------
% Time evolution 
% ----------------------------
for p = p
     [RX,RT,RXx, R2T,mtmTM, X10aux,X20aux, XX, TT, mtmT_q,mtmT,cpTaux,cpTaux2]= VFE_alg_ver3(M,MM,l,N,p,q) ;
     TTT = [TT(1, :); .5 * (TT(1:end-1, :) + TT(2:end, :)); TT(end, :)] ; 
     T = [TTT(:, 1)  TTT(:, 2) TTT(:,3)] ./ sqrt(TTT(:, 1).^2 - TTT(:, 2).^2- TTT(:, 3).^2) ; 
     X10(p+1,:) = X10aux ;   
     
    cpT(p+1,:) = cpTaux ;     
    cpT2(p+1,:) = cpTaux2 ;
    
    mtmTT(p+1,:) = mtmT ;
    mtmTTq(p+1,:) = mtmT_q ;
    mtmTTM(p+1,:) = mtmTM ;
%     
end

myfile = ['VFE_alg_M8_l06_Irr_p' num2str(p) 'q' num2str(q) '.mat'];
save(myfile, 'RT', 'RX', 'l','M','p','q');


toc 
    

