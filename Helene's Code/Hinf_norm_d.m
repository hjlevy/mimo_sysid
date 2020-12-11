function [H_infn,w_max] = Hinf_norm_d(A,B,C,D,gam,tol,ts)
%finding the Hinf norm and wmax of a discrete time system

%converting discrete time matrices to continuous
I = eye(size(A));
Ac = -inv(I+A)*(I-A);
Bc = sqrt(2)*inv(I+A)*B;
Cc = sqrt(2)*C*inv(I+A);
Dc = D - C*inv(I+A)*B;

[H_infn,nu] = Hinf_norm_c(Ac,Bc,Cc,Dc,gam,tol);

w_max = log((1 + 1j*nu)/(1 - 1j*nu))/(1j*ts);
w_max = real(w_max);
end