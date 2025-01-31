% -------------------------------------------------------
% The following function constructs the initial data i.e. the hyperbolic
% planar polygon
% Arguments: 
% M: Total number of sides
% N: The discretized nodes
% l = L /(M-1)
% -------------------------------------------------------

function [X,T] = GetPol(N,M,l,b)

% ---------------------------
% The tangent vectors
% ---------------------------
    a = sqrt(1+b^2) ;
    k = 0 : M-1 ; k=k.' ;
    T_vert = [a*cosh((k+(1-M)/2)*l) a*sinh((k+(1-M)/2)*l) b*ones(M,1)];
    
    for k = 0 : M-1
        T(k * N/M + (1:N/M),:) = ones(N/M,1) * T_vert(k+1,:) ; 
    end 

% -----------------------------------------------------
% drawing the polygon i.e. by integrating X_s = T 
% -----------------------------------------------------  

% with initial condition X(0)  
    
%     c20 = -l * M/2  ;
%     
%     X100 = l * 0.5 * sinh(c20 ) / sinh(l/2) ; 
%     X200 = l * 0.5 * cosh(c20 ) / sinh(l/2) ; 
%    
    L = M*l ; h = L / N ; 
%     T10 = T(:,1) ; T20 = T(:,2) ;
%     
%     X1 = Int_X(X100,T10,h) ; 
%     X2 = Int_X(X200,T20,h) ; 
%     
%     X = [X1  X2 zeros(N+1,1)] ; 
    X = h*cumsum(T); X = [0 0 0 ; X] ; 
    X1mean = h*0.5* (X(1,1) + 2*sum(X(2:end-1,1)) + X(end,1))/L ;
    X2mean = h*0.5* (X(1,2) + 2*sum(X(2:end-1,2)) + X(end,2))/L ;
    X3mean = h*0.5* (X(1,3) + 2*sum(X(2:end-1,3)) + X(end,3))/L ;
    X = X- [X1mean X2mean X3mean] ; 
%    Making T a N+1 point vector by taking the average of points

    TT = [T(1, :); .5 * (T(1:end-1, :) + T(2:end, :)); T(end, :)] ; 

    T = [TT(:, 1)  TT(:, 2) TT(:,3)] ./ sqrt(TT(:, 1).^2 - TT(:, 2).^2- TT(:, 3).^2) ; 
    
    
end


