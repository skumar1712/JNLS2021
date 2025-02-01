% ------------------------------------------------
% Symbolic calculation for rotation matrix 
% r = \rho, t = \theta
% Date: July 13, 2017
% The first part of the code evaluates the Matrix H = H_{q-1} * H_{q-2} ...
% * H_1 * H_0 by evaluating theta from the Gauss sums 
% The second part calculates the Rotation matrix by using the algorithm
% used in the paper 
% Minkowskii inner product : <x,y> = -x_1*y_1 + x_2*y_2 + x_3*y_3
% Hyperbolic case:
% Norm: || T || = -1,  || e_1 || = 1,  || e_2 || = 1

% Last update: August 7, 2019 (optimized a bit by removing for loops)
% ------------------------------------------------
% ------------------------------------------------
% close all; 
% clear 
% tic 
function [RX,RT,RXx, RTt,mtmTM, X10,X20, XX, TT,mtmT_q,mtmT,cpT,cpT2] = VFE_alg_ver3(M,MM,l,N,p1,q1) 

if mod(M,2)==1
    disp('Please choose M even!')
    error('Error!') 
end
L = l * M ; 
if p1 == 0 
    q = 1 ;
    p = 1 ; 
elseif ne(gcd(p1,q1),1) 
    p = p1 / gcd(p1,q1) ; 
    q = q1 / gcd(p1,q1) ; 
else 
    p = p1 ; 
    q = q1 ;
end
% ------------------------------------------------
% Computation of \theta_m using values of Gauss sums
% ------------------------------------------------
t = zeros(1,q) ; 
rhoq = zeros(1,q) ; 
if mod(q,2) == 1 
    for m = 0 : (q-1)/2
        t(m+1) = -1i*log(Gauss_sum(-p,m,q) / sqrt(q)) ; 
        rhoq(m+1) = acosh(2*cosh(l/2)^(2/q)-1) ;
    end
%     Using symmetries of Gauss sums
    t((q+3)/2:q) = fliplr(t(2:(q+1)/2)); 
    rhoq((q+3)/2:q) = fliplr(rhoq(2:(q+1)/2)); 
else
    for m = 0 : q/2
        if mod(q/2-m,2) == 0 
            t(m+1) = -1i*log(Gauss_sum(-p,m,q) / sqrt(2*q)) ; 
            rhoq(m+1) = acosh(2*cosh(l/2)^(4/q)-1) ;
        end
    end
%     Using symmetries of Gauss sums
    t(q/2+2:q) = fliplr(t(2:q/2)); 
    rhoq(q/2+2:q) = fliplr(rhoq(2:q/2)); 
end
t = real(t) ; 
% ------------------------------------------------
% Algebraic calculation for the hyperbolic case:
% ------------------------------------------------
% Matrix B stands for the Hyperbolic case
 if mod(q,2) == 1
    for j = 1 : (q+1)/2
        B(:,:,j)  = [cosh(rhoq(j)) sinh(rhoq(j))*cos(t(j)) sinh(rhoq(j))*sin(t(j)) ; ...
                sinh(rhoq(j))*cos(t(j)) cosh(rhoq(j))*cos(t(j))^2 + sin(t(j))^2 (cosh(rhoq(j))-1)*cos(t(j))*sin(t(j)); ...
                sinh(rhoq(j))*sin(t(j)) (cosh(rhoq(j))-1)*cos(t(j))*sin(t(j)) cosh(rhoq(j))*sin(t(j))^2 + cos(t(j))^2] ;
    end 
    B(:,:,(q+3)/2:q) = B(:,:,(q+1)/2:-1:2) ;
else
    for j = 1 : q/2+1
        B(:,:,j)  = [cosh(rhoq(j)) sinh(rhoq(j))*cos(t(j)) sinh(rhoq(j))*sin(t(j)) ; ...
                sinh(rhoq(j))*cos(t(j)) cosh(rhoq(j))*cos(t(j))^2 + sin(t(j))^2 (cosh(rhoq(j))-1)*cos(t(j))*sin(t(j)); ...
                sinh(rhoq(j))*sin(t(j)) (cosh(rhoq(j))-1)*cos(t(j))*sin(t(j)) cosh(rhoq(j))*sin(t(j))^2 + cos(t(j))^2] ;
    end
    B(:,:,q/2+2:q) = B(:,:,q/2:-1:2) ;
end
% Computation of Frenet frame {T,e1,e2} and tangent vectors 
Te1e2 = eye(3) ; 
T(1,:) = Te1e2(1,:) ; 
for j = 1 : M*q
    k = mod(j,q)  ;
    if k == 0 
        k = q ; 
    end
    Te1e2 = B(:,:,k) * Te1e2 ;
    T(j+1,:) = Te1e2(1,:) ; 
end 
T = T(2:end,:) ; % See eq. (45) and (52) T(1,:) = T(0+) which is obtained after the multiplication of R_1
% ------------------------------------------------
% Integrate X_s = T 

[R,X] = Get_X(T,M,l,q) ;
h = l/q ; 
RX = (R * X.').' ; 
RT = (R * T.').' ;
% RX(:,3) = RX(:,3) - mean(RX(:,3)) ; 
% RX(:,2) = RX(:,2) - mean(RX(:,2)) ; 
% RX(:,1) = RX(:,1) - mean(RX(:,1)) ; 
RX1mean = 0.5*(h/L)*(RX(1,1)+2*sum(RX(2:end-1,1))+RX(end,1)) ;
RX2mean = 0.5*(h/L)*(RX(1,2)+2*sum(RX(2:end-1,2))+RX(end,2)) ;
RX3mean = 0.5*(h/L)*(RX(1,3)+2*sum(RX(2:end-1,3))+RX(end,3)) ;
RX = RX - [RX1mean RX2mean RX3mean];

% RX = RX - mean(RX) ; 


% When MM/2 even 

for j = -(MM/2-1) : MM/2
    jj = j+(MM/2-1);
    RTt(jj*M*q/2+(1:M*q/2),:) = [cosh(j*l)*RT(1:M*q/2,1)+sinh(j*l)*RT(1:M*q/2,2) ...
                    sinh(j*l)*RT(1:M*q/2,1)+cosh(j*l)*RT(1:M*q/2,2) RT(1:M*q/2,3)] ;
end

[R,Xx] = Get_X(RTt,MM,l,q) ;

RXx = (R * Xx.').' ; 

LL = MM*l ;
RX1mean = 0.5*(h/LL)*(RXx(1,1)+2*sum(RXx(2:end-1,1))+RXx(end,1)) ;
RX2mean = 0.5*(h/LL)*(RXx(1,2)+2*sum(RXx(2:end-1,2))+RXx(end,2)) ;
RX3mean = 0.5*(h/LL)*(RXx(1,3)+2*sum(RXx(2:end-1,3))+RXx(end,3)) ;
RXx = RXx - [RX1mean RX2mean RX3mean];

sv = L*(0:M*q)/(M*q) - L/2 ; 

sMq = L*(0:M*q)/(M*q) - L/2 ;
i1x = floor((M*q)*(0:N-1)/N)+1 ;
i1t = floor((M*q)*(0:N-1)/N)+1 ;
s = L*(0:N)/N - L/2 ; 
s1 = sMq(i1x) ; 
i2x = floor((M*q)*(0:N-1)/N)+2 ;
s2 = sMq(i2x) ; 
XX = zeros(N+1,3); TT = zeros(N,3); 
for m = 1 : N/2+1
    XX(m,:) = (RX(i1x(m),:)*(s(m)-s2(m)) + RX(i2x(m),:)*(s1(m)-s(m))) / (s1(m)-s2(m)) ;    
end
XX(N/2+2:N+1,:) = [-flipud(XX(1:N/2,1)) flipud(XX(1:N/2,2)) flipud(XX(1:N/2,3))] ; 
for m = 1 : N
    TT(m,:) = RT(i1t(m),:); 
end

XX1mean = 0.5*(1/N)*(XX(1,1)+2*sum(XX(2:end-1,1))+XX(end,1)) ;
XX2mean = 0.5*(1/N)*(XX(1,2)+2*sum(XX(2:end-1,2))+XX(end,2)) ;
XX3mean = 0.5*(1/N)*(XX(1,3)+2*sum(XX(2:end-1,3))+XX(end,3)) ;
XX = XX - [XX1mean XX2mean XX3mean];
X10 = XX(N/2+1,:) ;
X20 = XX(N/2+1+N/(2*M),:) ;

mtmT = (l/q)*sum(cross(RX(1:M*q,:),RT(1:M*q,:)),1) ; mtmT(:,1) = -mtmT(:,1);
mtmTM = (l/q)*sum(cross(RXx(1:MM*q,:),RTt(1:MM*q,:)),1) ; mtmTM(:,1) = -mtmTM(:,1);

mtmT_q = (l/q)*sum(cross(RX(1:q,:),RT(1:q,:)),1); mtmT_q(:,1) = -mtmT_q(:,1);

cpT = sum(cross(RT(1:end-1,:),RT(2:end,:)),1)  ;
cpT2 = sum(cross(RTt(1:end-1,:),RTt(2:end,:)),1)  ;

end
