function [D,F] = FDmat(N,h)

%     Five points stencil
%     Both D and F are N x N matrices

    D11_ = [-25 48 -36 16 -3 zeros(1,N-5)]; 
    D1_ = [-3 -10 18 -6 1 zeros(1,N-5)] ; 
    D2_ = [1 -8 0 8 -1 zeros(1,N-5)]; 

    D1 = sparse(1:N-6,2:N-5, 1,N-6,N) ;
    D2 = sparse(1:N-6,3:N-4,-8,N-6,N) ;
%   D3 = sparse(1:N-6,4:N-3,0,N-6,N) ; 
    D4 = sparse(1:N-6,5:N-2,8,N-6,N) ; 
    D5 = sparse(1:N-6,6:N-1,-1,N-6,N) ; 
    D = D1 + D2 + D4 + D5; 
    Dn_ = [zeros(1,N-5) 1 -8 0 8 -1]; 
    Dn1_ = [zeros(1,N-5) -1 6 -18 10 3]; 
    Dnn1_ = [zeros(1,N-5) 3 -16 36 -48 25]; 
    D = [D11_; D1_ ; D2_ ; D ; Dn_ ; Dn1_; Dnn1_] /(12 * h) ;

    F11_ = [45 -154 214 -156 61 -10 zeros(1,N-6)] ; 
    F1_ = [10 -15 -4 14 -6 1 zeros(1,N-6)] ; 
    F2_ = [-1 16 -30 16 -1 zeros(1,N-5)] ; 

    F1 = sparse(1:N-6,2:N-5, -1,N-6,N) ;
    F2 = sparse(1:N-6,3:N-4,16,N-6,N) ;
    F3 = sparse(1:N-6,4:N-3,-30,N-6,N) ; 
    F4 = sparse(1:N-6,5:N-2,16,N-6,N) ; 
    F5 = sparse(1:N-6,6:N-1,-1,N-6,N) ; 
    F = F1 + F2 + F3 + F4 + F5; 
    
    Fn_ = [zeros(1,N-5) -1 16 -30 16 -1]; 
    Fn1_ = [zeros(1,N-6) 1 -6 14 -4 -15 10]; 
    Fn11_ = [zeros(1,N-6) -10 61 -156 214 -154 45]; 


    F = [F11_; F1_ ; F2_ ; F ; Fn_ ; Fn1_; Fn11_] /(12 * h^2) ;

return; 

% %     Seven points stencil 
%     D1_ = [-147 360 -450 400 -225 72 -10 zeros(1,N-7)] ; 
%     D2_ = [-10 -77 150 -100 50 -15 2 zeros(1,N-7)]; 
%     D3_ = [2 -24 -35 80 -30 8 -1 zeros(1,N-7)]; 
%     D4_ = [-1 9 -45 0 45 -9 1 zeros(1,N-7)] ;  
%     
%     D1 = sparse(1:N-8,2:N-7,-1,N-8,N) ;
%     D2 = sparse(1:N-8,3:N-6,9,N-8,N) ;
%     D3 = sparse(1:N-8,4:N-5,-45,N-8,N) ; 
%     
%     D4 = sparse(1:N-8,6:N-3,45,N-8,N) ; 
%     D5 = sparse(1:N-8,7:N-2,-9,N-8,N) ; 
%     D6 = sparse(1:N-8,8:N-1,1,N-8,N) ; 
%     D = D1 + D2 + D3 + D4 + D5 + D6; 
%     
%     Dn3_ = [zeros(1,N-7) -1 9 -45 0 45 -9 1] ;  
%     Dn2_ = [zeros(1,N-7) 1 -8 30 -80 35 24 -2]; 
%     Dn1_ = [zeros(1,N-7) -2 15 -50 100 -150 77 10];
%     Dn_ = [zeros(1,N-7) 10 -72 225 -400 450 -360 147] ; 
%     
%     D = [D1_ ; D2_ ;D3_; D4_; D ; Dn3_; Dn2_ ; Dn1_; Dn_] /(60 * h) ;
%     
%     F1_ = [812 -3132 5265 -5080 2970 -972 137 zeros(1,N-7)] ; 
%     F2_ = [126  -70  -486  855  -670  324  -90  11 zeros(1,N-8)] ; % 8 point scheme
% %     F2_ = [137 -147 -255 470 -285 93 -13 zeros(1,N-7)]; 
% %     F3_ = [-13 228 -420 200 15 -12 2 zeros(1,N-7)]; 
%     F3_ = [-11 214 -378 130  85 -54 16 -2 zeros(1,N-8)] ; 
%     F4_ = [2 -27 270 -490 270 -27 2 zeros(1,N-7)] ;  
%     
%     F1 = sparse(1:N-8,2:N-7,2,N-8,N) ;
%     F2 = sparse(1:N-8,3:N-6,-27,N-8,N) ;
%     F3 = sparse(1:N-8,4:N-5,270,N-8,N) ; 
%     F4 = sparse(1:N-8,5:N-4,-490,N-8,N) ; 
%     F5 = sparse(1:N-8,6:N-3,270,N-8,N) ; 
%     F6 = sparse(1:N-8,7:N-2,-27,N-8,N) ; 
%     F7 = sparse(1:N-8,8:N-1,2,N-8,N) ; 
%     F = F1 + F2 + F3 + F4 + F5 + F6 + F7; 
%     
%     Fn3_ = [zeros(1,N-7) 2 -27 270 -490 270 -27 2] ;  
% %     Fn2_ = [zeros(1,N-7) 2 -12 15 200 -420 228 -13]; 
%     Fn2_ = fliplr(F3_) ; 
% %     Fn1_ = [zeros(1,N-7) -13 93 -285 470 -255 -147 137];
%     Fn1_ = fliplr(F2_) ; 
%     Fn_ = [zeros(1,N-7) 137 -972 2970 -5080 5265 -3132 812] ; 
%     
%     F = [F1_ ; F2_ ;F3_; F4_; F ; Fn3_; Fn2_ ; Fn1_; Fn_] /(180 * h^2) ;
%     return; 
    
% 
% 
%     F1_ = [11 -20 6 4 -1 zeros(1,N-5)] ; 
%     F2_ = [-1 16 -30 16 -1 zeros(1,N-5)] ; 
% 
%     F1 = sparse(1:N-6,2:N-5, -1,N-6,N) ;
%     F2 = sparse(1:N-6,3:N-4,16,N-6,N) ;
%     F3 = sparse(1:N-6,4:N-3,-30,N-6,N) ; 
%     F4 = sparse(1:N-6,5:N-2,16,N-6,N) ; 
%     F5 = sparse(1:N-6,6:N-1,-1,N-6,N) ; 
%     F = F1 + F2 + F3 + F4 + F5; 
%     Fn_ = [zeros(1,N-5) -1 16 -30 16 -1]; 
%     Fn1_ = [zeros(1,N-5) -1 4 6 -20 11]; 
%     F = [F1_ ; F2_ ; F ; Fn_ ; Fn1_] /(12 * h^2) ;
%     
    
%     return ; 

    
    
    
    
    
    % D1_ = [-25 48 -36 16 -3 zeros(1,N-5)] ; 
% D2_ = [-3 -10 18 -6 1 zeros(1,N-5)]; 
% 
% D1 = sparse(1:N-4,1:N-4, 1,N-4,N) ;
% D2 = sparse(1:N-4,2:N-3,-8,N-4,N) ;
% D3 = sparse(1:N-4,3:N-2,0,N-4,N) ; 
% D4 = sparse(1:N-4,4:N-1,8,N-4,N) ; 
% D5 = sparse(1:N-4,5:N,-1,N-4,N) ; 
% D = D1 + D2 + D3 + D4 + D5; 
% Dn_ = [zeros(1,N-5) -1 6 -18 10 3]; 
% Dn1_ = [zeros(1,N-5) 3 -16 36 -48 25]; 
% D = [D1_ ; D2_ ; D ; Dn_ ; Dn1_] /(12 * h) ;


% F = d^2/dx^2



% F1_ = [35 -104 114 -56 11 zeros(1,N-5)] ; 
% F2_ = [11 -20 6 4 -1 zeros(1,N-5)]; 
% 
% F1 = sparse(1:N-4,1:N-4, -1,N-4,N) ;
% F2 = sparse(1:N-4,2:N-3,16,N-4,N) ;
% F3 = sparse(1:N-4,3:N-2,-30,N-4,N) ; 
% F4 = sparse(1:N-4,4:N-1,16,N-4,N) ; 
% F5 = sparse(1:N-4,5:N,-1,N-4,N) ; 
% F = F1 + F2 + F3 + F4 + F5; 
% Fn_ = [zeros(1,N-5) -1 4 6 -20 11]; 
% Fn1_ = [zeros(1,N-5) 11 -56 114 -104 35]; 
% F = [F1_ ; F2_ ; F ; Fn_ ; Fn1_] /(12 * h^2) ;

