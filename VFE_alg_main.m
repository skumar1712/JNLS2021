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
% q=2^3*3^2*5^2*7; 
% for jj = 5
% q = 2^(jj-1)*500+1  ;
M= 8; 
p = 10327 ; q = 27993 ;% for p/q = 1/3 + 1/31 + 1/301  
p = 18209 ; q = 65764 ; 
p = 6004; q = 18003; % p/q = 1/3 + 1/3001
% p = 5004 ; q = 15003; % p/q = 1/3 + 1/5001
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
%      RX : corresponds to M
%      RXx: corresponds to MM 

%      [RX,RT,X10aux, XX, TT,mtmX,mtmX_q,mtmT_q,mtmT,cpT] = VFE_alg_ver2(M,l,N,p,q) ;
     TTT = [TT(1, :); .5 * (TT(1:end-1, :) + TT(2:end, :)); TT(end, :)] ; 
     T = [TTT(:, 1)  TTT(:, 2) TTT(:,3)] ./ sqrt(TTT(:, 1).^2 - TTT(:, 2).^2- TTT(:, 3).^2) ; 
%      X1fu_alg(:,p+1) = XX(:,1); X2fu_alg(:,p+1) = XX(:,2); X3fu_alg(:,p+1) = XX(:,3); 
%      T1fu_alg(:,p+1) = T(:,1); T2fu_alg(:,p+1) = T(:,2); T3fu_alg(:,p+1) = T(:,3); 
     X10(p+1,:) = X10aux ;   
     
%      T00(p+1,:) = T(N/2,:) ; 
%      X20(p+1,:) = X20aux ;    
%      
     
%     X1fu_alg(:,p+1) = RX(:,1); X2fu_alg(:,p+1) = RX(:,2); X3fu_alg(:,p+1) = RX(:,3); 
%     T1fu_alg(:,p+1) = RT(:,1); T2fu_alg(:,p+1) = RT(:,2); T3fu_alg(:,p+1) = RT(:,3); 
    cpT(p+1,:) = cpTaux ;     
    cpT2(p+1,:) = cpTaux2 ;
    
    mtmTT(p+1,:) = mtmT ;
    mtmTTq(p+1,:) = mtmT_q ;
    mtmTTM(p+1,:) = mtmTM ;
%     
end

myfile = ['VFE_alg_M8_l06_Irr_p' num2str(p) 'q' num2str(q) '.mat'];
save(myfile, 'RT', 'RX', 'l','M','p','q');

% 
% ang=  l/2; 
% 
% R = [cosh(ang) -sinh(ang) ; -sinh(ang) cosh(ang)] ; 
% X20e = X20-X20(end,:);
% X20er = [(R*X20e(:,1:2).').' X20e(:,3)];
% X10a = X10 - X10(1,:) ;
% z1 = X10a(:,2)+1i*X10a(:,3) ;
% z2 = X20er(:,2)+1i*X20er(:,3) ;
% z3= (conj(z1(end))+z2) ; 
% z4 = (conj(conj(z1(end))+z2)) ; 
% ReZ = [real(z1) ; real(z3(end:-1:1)); real(z4); real(z1(end:-1:1)) ];
% ImZ = [imag(z1) ; -imag(z3(end:-1:1)); -imag(z4); -imag(z1(end:-1:1)) ];
% 
%     

% savefile = ['VFE_alg_M' num2str(M) '_2N' num2str(log2(N/M)) 'l00' num2str(l*1e3) 'q' num2str(q) '.mat'];


% save(savefile, 'M', 'l', 'X1fu_alg', 'X2fu_alg', 'X3fu_alg','T1fu_alg', 'T2fu_alg', 'T3fu_alg','X0','-v7.3') 
% save(savefile, 'M', 'l', 'mtmXX', 'mtmXXq', 'mtmTT','mtmTTq','cpT','X0','-v7.3') 
% save(savefile, 'M', 'l', 'RT', 'RX', '-v7.3') 
% end

toc 
    

