% ---------------------------------------------------------------
% The program solves VFE: X_t = X_s \wedege X_{ss} only for X(s,0)
% a hyperbolic polygon parameterized by a parameter l > 0 and number of
% sides given by M
% Date: February 01, 2019

% Last modified: July 30, 2019 
% ---------------------------------------------------------------

clear
tic 


n = 7 ;
l_arr = [ 0.4 ] ; 

    for jj = 1 : length(l_arr) 
        
% Parameters

% n = 5 ; 

M = 48;
N =  2^n * M ; 
% l = 0.025 ; 
l = l_arr(jj) ;
T_run = l^2/(2*pi); 
c0 = sqrt(2*log(cosh(l/2))/pi) ; 
L = l * M ; 
h = L / N ; 
b = 0.4 ;   % any number 0<=b<inf
a = sqrt(1+b^2) ;
q = 1260 ; 
m = 0 ;
% N_t = 4^(n-7)*4980;
N_t = 4^(n-7)*5040;
N_t = 2^m*round(N_t) ; 
dt = T_run / N_t ; 
t = linspace(0,T_run,N_t+1);

% ---------------------------------------------------------------
% Stability condition
if dt / h^2 > 0.5298 % (best) 
    error('Step size is too large!')
end
% -----------------------------------------------------     
% ----------------MAT FILE NAME 
savefile = ['VFE_BCT_2N' num2str(n) 'M' num2str(M) 'l00' num2str(l*1e3) 'q' num2str(q) '.mat']; 
% savefile = sprintf('VFE_BCT_2N%dM%d_l0025.mat',n,M)  ;

% ---------------------------------------------------------------
% Differentiation matrices 
% d/dx and d^2/dx^2 matrices

[DD,FF] = FDmat(N+1,h) ; % D, F are (N+1)x(N+1) dim
D = DD(2:end-1,:) ; 
F = FF(2:end-1,:) ; 
s = linspace(-L/2,L/2,N+1) ; 

% ---------------------------------------------------------------
% Initial data for X and T

[X,T] = GetPol(N,M,l,b) ; 
% T = [a*cosh(s).' a*sinh(s).' b*ones(N+1,1)];
% X = [a*sinh(s).' a*cosh(s).' b*s.']  ; 

X10 = X(:,1) ; X20 = X(:,2) ;  X30 = X(:,3) ; 
X1 = X10 ; X2 = X20 ; X3 = X30; 
T10 = T(:,1) ; T20 = T(:,2); T30 = T(:,3) ; 
T1 = T10 ; T2 = T20 ; T3 = T30; 

% ---------------------------------------------------------------
% Preallocation 

X1mean = zeros(1,N_t+1) ; X2mean = zeros(N_t+1,M/4+1) ; 
X3mean = zeros(N_t+1,M/4+1) ; X2meana = zeros(N_t+1,M/4+1) ; 
% X2mean_N = zeros(1,N_t+1) ; X3mean_N = zeros(1,N_t+1) ; 
% X2mean_inn = zeros(1,N_t+1) ; X2mean_s1 = zeros(1,N_t+1); X2mean_s2 = zeros(1,N_t+1); 
% X3mean_inn = zeros(1,N_t+1) ; X3mean_s1 = zeros(1,N_t+1); X3mean_s2 = zeros(1,N_t+1); 
% X2mean_1(1) = mean(X2(1:N/2)); 
% X2mean_2(1) = mean(X2(N/2+2:end)); 
% X2mean_3(1) = mean(X2(1:end-1)); 
% X2mean_inn1(1) = mean(X2(N/4+1:N/2)); 

X2tmean = zeros(1,N_t+1) ; X3tmean = zeros(1,N_t+1) ; 

z = zeros(N_t+1,3) ; 
X1full = zeros(N+1,q); X2full = zeros(N+1,q); X3full = zeros(N+1,q); 
T1full = zeros(N+1,q); T2full = zeros(N+1,q); T3full = zeros(N+1,q); 
T_norm = zeros(N_t+1,1) ; 

% ---------------------------------------------------------------
% Storing the initial data 

X1full(:,1) = X1 ; X2full(:,1) = X2 ; X3full(:,1) = X3 ; 
T1full(:,1) = T1 ; T2full(:,1) = T2 ; T3full(:,1) = T3 ; 
Tnorm = max(sqrt(T1.^2 - T2.^2 - T3.^2)) ; 
T_norm(1) = Tnorm; 

X1mean(1) = mean(X1) ; 

%     Center of mass 
    for r = 0 : M/4-1
        X3mean(1,r+1) = mean(X3(2*r*N/M+1:end-(2*r*N/M+1))); 
        X2mean(1,r+1) = mean(X2(2*r*N/M+1:end-(2*r*N/M+1)));
        X2meana(1,r+1) = mean(X2(2*r*N/M+1:end-(2*r*N/M)));
    end
    X3mean(1,M/4+1) = mean(X3(N/2+1-N/M:N/2+N/M)); 
    X2mean(1,M/4+1) = mean(X2(N/2+1-N/M:N/2+N/M)); 
    X2meana(1,M/4+1) = mean(X2(N/2+1-N/M:N/2+N/M+1)); 
% 
%     X3mean(j+1,1) = mean(X3(N/2+1+(1:N/M))); 
%     X2mean(j+1,1) = mean(X2(N/2+1+(1:N/M))); 
%     for r = 1 : M/4
%         X3mean(j+1,r+1) = mean(X3(N/2+1+(1:2*r*N/M))); 
%         X2mean(j+1,r+1) = mean(X2(N/2+1+(1:2*r*N/M))); 
%     end
%     Third component:
%     X3mean_N(1) = mean(X3) ;
%     X3mean_1(j+1) = mean(X3(1:N/2)); 
%     X3mean_inn(j+1) = mean(X3(N/4+1:end-N/4)) ;
%     X3mean_s1(j+1) = mean(X3(N/2+1:N/2+N/M)) ;  % one side
%     X3mean_s2(j+1) = mean(X3(N/4+1:N/4+N/M)) ;  % one side
%     
% %    Second component:
%     X2mean_N(1) = mean(X2) ;
%     X2mean_inn(j+1) = mean(X2(N/4+1:end-N/4)) ;
%     X2mean_s1(j+1) = mean(X2(N/2+1:N/2+N/M)) ;  % one side
%     X2mean_s2(j+1) = mean(X2(N/4+1:N/4+N/M)) ;  % one side

z(1,:) = [X1(N/2+1) X2(N/2+1) X3(N/2+1)] ;

% -----------------------------------------------------------------
    T1s = DD * T1 ; T2s = DD * T2 ; T3s = DD * T3 ;
    AX1 = -(T2 .* T3s - T3 .* T2s) ; AX2 = -(T1 .* T3s - T3 .* T1s) ; 
    AX3 = (T1 .* T2s - T2 .* T1s) ; 

    X2tmean(1) = mean(AX2(N/4+1:end-N/4)) ;
    X3tmean(1) = mean(AX3(N/4+1:end-N/4))  ; 

% ---------------------------------------------------------------
% Time evolution

p=1; 

for j = 1 : N_t 
    
    if max(abs(T1))>cosh(L)
    	j
        error('Error! Norm exceeded!') 
    end 
    
    if isnan(X1)
        error('Error!')
    end
    
    T1ss = F * T1; T1s = DD * T1 ;
    T2ss = F * T2; T2s = DD * T2 ;
    T3ss = F * T3; T3s = DD * T3 ;
    
    AT1 = -(T2(2:end-1) .* T3ss - T3(2:end-1) .* T2ss) ; 
    AT2 = -(T1(2:end-1) .* T3ss - T3(2:end-1) .* T1ss) ;
    AT3 = (T1(2:end-1) .* T2ss - T2(2:end-1) .* T1ss) ; 
    
    AX1 = -(T2 .* T3s - T3 .* T2s) ; 
    AX2 = -(T1 .* T3s - T3 .* T1s) ;
    AX3 = (T1 .* T2s - T2 .* T1s) ; 
    
    T1_aux = T1 + 0.5 * dt * [0 ; AT1; 0] ; 
    T2_aux = T2 + 0.5 * dt * [0 ; AT2; 0] ; 
    T3_aux = T3 + 0.5 * dt * [0 ; AT3; 0] ;  
    Tnorm = sqrt(T1_aux.^2-T2_aux.^2-T3_aux.^2) ; 
    T1_aux = T1_aux ./ Tnorm ; T2_aux = T2_aux ./ Tnorm;  T3_aux = T3_aux ./ Tnorm; 
    
    T1ss = F * T1_aux ; T1s = DD * T1_aux ;
    T2ss = F * T2_aux ; T2s = DD * T2_aux ;
    T3ss = F * T3_aux ; T3s = DD * T3_aux ;
    
    BT1 = -(T2_aux(2:end-1) .* T3ss - T3_aux(2:end-1) .* T2ss) ; 
    BT2 = -(T1_aux(2:end-1) .* T3ss - T3_aux(2:end-1) .* T1ss) ;
    BT3 = (T1_aux(2:end-1) .* T2ss - T2_aux(2:end-1) .* T1ss) ; 
    
    BX1 = -(T2_aux .* T3s - T3_aux .* T2s) ; 
    BX2 = -(T1_aux .* T3s - T3_aux .* T1s) ;
    BX3 = (T1_aux .* T2s - T2_aux .* T1s) ; 
    
    T1_aux = T1 + 0.5 * dt * [0 ; BT1; 0] ; 
    T2_aux = T2 + 0.5 * dt * [0 ; BT2; 0] ; 
    T3_aux = T3 + 0.5 * dt * [0 ; BT3; 0] ;  
    Tnorm = sqrt(T1_aux.^2-T2_aux.^2-T3_aux.^2) ; 
    T1_aux = T1_aux ./ Tnorm ; T2_aux = T2_aux ./ Tnorm;  T3_aux = T3_aux ./ Tnorm; 
    
    T1ss = F * T1_aux ; T1s = DD * T1_aux ;
    T2ss = F * T2_aux ; T2s = DD * T2_aux ;
    T3ss = F * T3_aux ; T3s = DD * T3_aux ;
    
    CT1 = -(T2_aux(2:end-1) .* T3ss - T3_aux(2:end-1) .* T2ss) ; 
    CT2 = -(T1_aux(2:end-1) .* T3ss - T3_aux(2:end-1) .* T1ss) ;
    CT3 = (T1_aux(2:end-1) .* T2ss - T2_aux(2:end-1) .* T1ss) ; 
    
    CX1 = -(T2_aux .* T3s - T3_aux .* T2s) ; 
    CX2 = -(T1_aux .* T3s - T3_aux .* T1s) ;
    CX3 = (T1_aux .* T2s - T2_aux .* T1s) ; 
    
    T1_aux = T1 +  dt * [0 ; CT1; 0] ; 
    T2_aux = T2 +  dt * [0 ; CT2; 0] ; 
    T3_aux = T3 +  dt * [0 ; CT3; 0] ;  
    Tnorm = sqrt(T1_aux.^2-T2_aux.^2-T3_aux.^2) ; 
    T1_aux = T1_aux ./ Tnorm ; T2_aux = T2_aux ./ Tnorm;  T3_aux = T3_aux ./ Tnorm; 
    
    T1ss = F * T1_aux ; T1s = DD * T1_aux ;
    T2ss = F * T2_aux ; T2s = DD * T2_aux ;
    T3ss = F * T3_aux ; T3s = DD * T3_aux ;
    
    DT1 = -(T2_aux(2:end-1) .* T3ss - T3_aux(2:end-1) .* T2ss) ; 
    DT2 = -(T1_aux(2:end-1) .* T3ss - T3_aux(2:end-1) .* T1ss) ;
    DT3 = (T1_aux(2:end-1) .* T2ss - T2_aux(2:end-1) .* T1ss) ; 
    
    DX1 = -(T2_aux .* T3s - T3_aux .* T2s) ; 
    DX2 = -(T1_aux .* T3s - T3_aux .* T1s) ;
    DX3 = (T1_aux .* T2s - T2_aux .* T1s) ; 
    
    RHSX1 = (AX1 + 2 * BX1 + 2 * CX1 + DX1)/ 6;
    RHSX2 = (AX2 + 2 * BX2 + 2 * CX2 + DX2)/6 ;
    RHSX3 = (AX3 + 2 * BX3 + 2 * CX3 + DX3)/6; 
    
%     X1tmean(j+1) = mean(RHSX1(4*N/M:end-4*N/M+1))  ; 
%     X2tmean(j+1) = mean(RHSX2(4*N/M:end-4*N/M+1)) ;
%     X3tmean(j+1) = mean(RHSX3(4*N/M:end-4*N/M+1))  ; 
    X2tmean(j+1) = mean(RHSX2(N/4+1:end-N/4)) ;
    X3tmean(j+1) = mean(RHSX3(N/4+1:end-N/4))  ; 
    
    X1 = X1 + dt *  (AX1 + 2 * BX1 + 2 * CX1 + DX1) / 6 ; 
    X2 = X2 + dt *  (AX2 + 2 * BX2 + 2 * CX2 + DX2) / 6  ; 
    X3 = X3 + dt *  (AX3 + 2 * BX3 + 2 * CX3 + DX3) / 6 ;
    
    T1 = T1 + dt * [0 ;(AT1 + 2 * BT1 + 2 * CT1 + DT1) / 6 ; 0] ;
    T2 = T2 + dt * [0 ;(AT2 + 2 * BT2 + 2 * CT2 + DT2) / 6 ; 0];
    T3 = T3 + dt * [0 ;(AT3 + 2 * BT3 + 2 * CT3 + DT3) / 6 ; 0];
    
    Tnorm = sqrt(T1.^2-T2.^2-T3.^2) ; 
    T1 = T1 ./ Tnorm ; T2 = T2 ./ Tnorm;  T3 = T3 ./ Tnorm; 
    
%     T_norm(j+1) = max(Tnorm) ; 

    if mod(j,N_t*p/q) == 0 
        X1full(:,p+1) = X1 ; X2full(:,p+1) = X2 ; X3full(:,p+1) = X3 ; 
        T1full(:,p+1) = T1 ; T2full(:,p+1) = T2 ; T3full(:,p+1) = T3 ; 
        p = p + 1  ;
        p/q ;
    end 
    
    z(j+1,:) = [X1(N/2+1) X2(N/2+1) X3(N/2+1)] ;
    X1mean(j+1) = mean(X1) ;
%     X1mean_inn(j+1) = mean(X1(N/4+1:end-N/4)) ;
%     Center of mass for inner points
%     X3mean(j+1,1) = mean(X3(N/2+1+(1:N/M))); 
    
    for r = 0 : M/4-1
        X3mean(j+1,r+1) = mean(X3(2*r*N/M+1:end-(2*r*N/M+1))); 
        X2mean(j+1,r+1) = mean(X2(2*r*N/M+1:end-(2*r*N/M+1)));
        X2meana(j+1,r+1) = mean(X2(2*r*N/M+1:end-(2*r*N/M)));
%         X2mean(j+1,r+1) = mean(X2(N/2+1+(1:2*r*N/M))); 
    end
    X3mean(j+1,M/4+1) = mean(X3(N/2+1-N/M:N/2+N/M)); 
%     X3mean(j+1,M/4+2) = mean(X3(N/2+1-N/M:N/2)); 
    X2mean(j+1,M/4+1) = mean(X2(N/2+1-N/M:N/2+N/M));
    X2meana(j+1,M/4+1) = mean(X2(N/2+1-N/M:N/2+N/M+1));
%     X2mean(j+1,M/4+2) = mean(X2(N/2+1-N/M:N/2)); 
    
%     X3mean_N(j+1) = mean(X3) ; 
%     X3mean_1(j+1) = mean(X3(1:N/2)); 
%     X3mean_2(j+1) = mean(X3(1:N/2+1)); 
%     X3mean_inn(j+1) = mean(X3(N/4+1:end-N/4)) ;
%     X3mean_inn1(j+1) = mean(X3(N/4+1:N/2)) ;
%     X3mean_inn2(j+1) = mean(X3(N/2+2:3*N/4+1)) ;
%     X3mean_s1(j+1) = mean(X3(N/2+1:N/2+N/M)) ; %    For one side
%     X3mean_s2(j+1) = mean(X3(N/4+1:N/4+N/M)) ; 
% %     X3_mean_s3(j+1) = mean(X3(1:N/M)) ;
%     
%     X2mean_N(j+1) = mean(X2) ; 
%     X2mean_inn(j+1) = mean(X2(N/4+1:end-N/4)) ;
%     
%     X2mean_1(j+1) = mean(X2(1:N/2)); 
%     X2mean_2(j+1) = mean(X2(N/2+2:end)); 
%     X2mean_3(j+1) = mean(X2(1:end-1)); 
%     X2mean_inn1(j+1) = mean(X2(N/4+1:N/2)) ;
%     X2mean_s1(j+1) = mean(X2(N/2+1:N/2+N/M)) ; 
%     X2mean_s2(j+1) = mean(X2(N/4+1:N/4+N/M)) ; 
%     ind = N/2-N/(2*M)+1: N/2+N/(2*M) ; 
%     X2_mean_b(j+1) = mean(X2(ind)) ;
    
    
end 

% c_M = (X3mean_N(end)-X3mean_N(1)) / T_run ;
% c_M_a = (X3mean_inn(end)-X3mean_inn(1)) / T_run ;
% c_M_s1 = (X3mean_s1(end)-X3mean_s1(1)) / T_run ;
% c_M_s2 = (X3mean_s2(end)-X3mean_s2(1)) / T_run ;
% c_M_s3 = (X3_mean_s3(end)-X3_mean_s3(1)) / T_run ;

% save(savefile,'X1','X2','X3','T1','T2','T3','dt','L') ; 

% error_h = norm(X3_mean - c_m * t,2) ; 
    c0 = sqrt(2*log(cosh(l/2))/pi);
    cM=2*pi*c0.^2 ./ (l .* sqrt(1-exp(-pi*c0.^2)));
    
% save(savefile, 'T1full',  'T2full', 'T3full', 'X1full', 'X2full', 'X3full', 'z',...
%      'X2mean_N', 'X2mean_inn','X2mean_s1','X2mean_s2', 'X3mean_N','X3mean_inn', 'X3mean_s1', 'X3mean_s2', ...
%      'X2tmean', 'X3tmean','l', 'dt', 'M', '-v7.3') 
save(savefile, 'T1full',  'T2full', 'T3full', 'X1full', 'X2full', 'X3full', 'z',...
     'X2mean','X2meana', 'X3mean','l', 'dt', 'M', '-v7.3')     

% X1_ex = a*sinh(s+b*t(end)).' ;
% X2_ex = a*cosh(s+b*t(end)).' ;
% X3_ex = b*s.'+a^2*t(end) ;
% 
% err = sqrt(mean((X1-X1_ex).^2 + (X2-X2_ex).^2 + (X3-X3_ex).^2)); 
    end
%     clear 
toc

return; 

% clear 
% ------------------------------------------------
% COMPUTING TANGENT VECTOR FROM N VALUES ------
% ------------------------------------------------

n = 8;
L = l*M; 
T1 = T1full(:,end) ; T2 = T2full(:,end) ; T3 = T3full(:,end) ; 
N = length(T3full)-1;
s = linspace(-L/2,L/2,N+1);
NM = N/M ;
for j = 1 : M/2
%     CASE 1: USING MEAN OF ONLY THE INNER N/2 POINTS
    ind = (NM/4+1 : 3*NM/4)+(j-1)*NM ; 
    T1v(j) = mean(T1(ind));
    T2v(j) = mean(T2(ind));
    T3v(j) = mean(T3(ind));    
%     CASE 2: USING MEAN OF ALL THE POINTS
%     T1v(j) = mean(T1((j-1)*2^n+1:j*2^n));
%     T2v(j) = mean(T2((j-1)*2^n+1:j*2^n));
%     T3v(j) = mean(T3((j-1)*2^n+1:j*2^n));
    sv(j) = mean(s((j-1)*2^n+1:j*2^n));
%     CASE 3: BY TAKING ONLY THE MIDDLE POINT 
% T3v(j) = T3(2^(n-1)+(j-1)*2^n+1) ;
% sv(j) = s(2^(n-1)+(j-1)*2^n+1) ;
end
for j = M/2+1:M
    ind = (NM/4+1 : 3*NM/4)+(j-1)*NM+1 ; 
    T1v(j) = mean(T1(ind));
    T2v(j) = mean(T2(ind));
    T3v(j) = mean(T3(ind));
%     T1v(j) = mean(T1((j-1)*2^n+2:j*2^n+1));
%     T2v(j) = mean(T2((j-1)*2^n+2:j*2^n+1));
%     T3v(j) = mean(T3((j-1)*2^n+2:j*2^n+1));
    sv(j) = mean(s((j-1)*2^n+2:j*2^n+1));
% T3v(j) = T3(2^(n-1)+(j-1)*2^n+1);
% sv(j) = s(2^(n-1)+(j-1)*2^n+1) ;
end
Tv = [T1v.' T2v.' T3v.'];
Tv=Tv./sqrt(Tv(:,1).^2-Tv(:,2).^2-Tv(:,3).^2);


toc
        
