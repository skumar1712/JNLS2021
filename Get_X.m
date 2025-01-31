function [R,X] = Get_X(T,M,l,q)

X(1,:) = [0 0 0] ; 
h = l/ (q) ; 
for j = 1 : M*q
    X(j+1,:) = X(j,:) + h * T(j,:) ; 
end
% ------------------------------------------------
% Rotation algorithm to align X in the correct way 
% ------------------------------------------------
% First rotation: to align into a plane orthogonal to z-axis 
% X = X(1:end-1,:) ;

ind1 = M*q/2 - q + 1 ; 
ind2 = M*q/2 + q + 1 ; 
vm = (X(ind1,:) - X(M*q/2+1,:)) ; 
vp = (X(ind2,:) - X(M*q/2+1,:)) ; 
% vp = X(end,:) - X(M_new/2+1,:) ; 
% vm = X(1,:) - X(M_new/2+1,:) ; 
% Normalisation since they are time-like 
vp = vp / sqrt(vp(1)^2-vp(2)^2-vp(3)^2) ;
vm = vm / sqrt(vm(1)^2-vm(2)^2-vm(3)^2) ;
v12=cross(vp,vm) ; v12(1) = -v12(1) ; % a space-like vector
% v12=cross(vm,vp) ; v12(1) = -v12(1) ; 
v12 = v12 / sqrt(-v12(1)^2+v12(2)^2+v12(3)^2) ; 
% z axis is a space-like vector so if v12 is space like vector 
% and if they span a time-like vector space i.e. 
% v12 . [0 0 1] > |v12| | 0 0 1| 
% then ... 
%  and the required axis of rotation is computed using the cross product of v12
% and [0 0 1]

if  v12(3)>1
    ang1 = acosh(v12(3)) ; 
    cp = cross(v12,[0 0 1]) ; cp(1) = -cp(1) ; 
    ax1 = cp / sqrt(-cp(1)^2+cp(2)^2+cp(3)^2) ;  % a space-like vector 
    K1 = [0 ax1(3) ax1(2) ; ax1(3) 0 ax1(1) ; ax1(2) -ax1(1) 0 ] ; % Semi-skew symmetric 
    R1 = eye(3) + sinh(ang1) * K1 + (cosh(ang1)-1) * K1^2 ; 
% else if they span space-like vector space i.e. 
% v12 . [0 0 1] < |v12| | 0 0 1| 
elseif  v12(3)<1
    if v12(3) == -1
        R1 = [-1 0 0; 0 -1 0; 0 0 1] ; 
    else
        ang1 = acos(v12(3)) ; 
        cp = cross(v12,[0 0 1]) ; cp(1) = -cp(1) ; 
        ax1 = cp / sqrt(cp(1)^2-cp(2)^2-cp(3)^2) ;  % a time-like vector 
        K1 = [0 ax1(3) ax1(2) ; ax1(3) 0 ax1(1) ; ax1(2) -ax1(1) 0 ] ; % Semi-skew symmetric 
        R1 = eye(3) + sin(ang1) * K1 + (1-cos(ang1)) * K1^2 ; 
    end
elseif  v12(3) == 1
%     ang1 = 0 ;     
    R1 = eye(3) ;     
end
% Second rotation to align vp - vm = c * (1,0,0) 
% R1X = (R1*X.').' ; 
% Type 1: 
R1vp = R1*vp.' ; R1vm = R1*vm.' ; 
% Type 2:
% R1vp = R1X(ind1,:) - R1X(M_new/2+1,:) ;  % s = -L/M
% R1vm = R1X(ind2,:) - R1X(M_new/2+1,:) ;  % s = L/M 
% Type 3: 
% % R1vp = R1X(end,:) - R1X(M_new/2+1,:) ; 
% % R1vm = R1X(1,:) - R1X(M_new/2+1,:) ;
vpm = R1vp - R1vm ;
vpm = vpm / sqrt(vpm(1)^2 - vpm(2)^2 - vpm(3)^2)  ;
if -vpm(1)^2 + vpm(2)^2 + vpm(3)^2 < 0 % time-like vector 
    ang2 = acosh(vpm(1));
    if abs(ang2) < 1e-10
        R2 = eye(3) ;
    else
        cp1 = cross(vpm,[1 0 0 ]); cp1(1) = -cp1(1) ; 
        ax2 = cp1 / sqrt(-cp1(1)^2+cp1(2)^2+cp1(3)^2 ) ; 
        K2 = [0 ax2(3) ax2(2) ; ax2(3) 0 ax2(1) ; ax2(2) -ax2(1) 0 ] ; % Semi-skew symmetric 
        R2 = eye(3) + sinh(ang2) * K2 + (cosh(ang2)-1) * K2^2 ; 
    end
% angle between space-like and time-like vector 
elseif -vpm(1)^2 + vpm(2)^2 + vpm(3)^2 > 0 
    ang2 = asinh(vpm(1)); 
    cp1 = cross(vpm,[1 0 0 ]); cp1(1) = -cp1(1) ; 
    if -cp1(1)^2+cp1(2)^2+cp1(3)^2 < 0 
        ax2 = cp1 / sqrt(cp1(1)^2-cp1(2)^2-cp1(3)^2 ) ; 
    else 
        ax2 = cp1 / sqrt(-cp1(1)^2+cp1(2)^2+cp1(3)^2 ) ; 
        K2 = [0 ax2(3) ax2(2) ; ax2(3) 0 ax2(1) ; ax2(2) -ax2(1) 0 ] ; % Semi-skew symmetric 
        R2 = eye(3) + sinh(ang2) * K2 + (cosh(ang2)-1) * K2^2 ; 
    end
else
    R2 = eye(3) ;
end
% Final Rotation matrix 
ang2
ax2
 R2
R = R2 * R1 ; 





