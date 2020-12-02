function [H_infn,w_max] = Hinf_norm_c(A,B,C,D,gam,tol)
Dg = gam^2*I-D'*D;
Aclp = [A+B*inv(Dg)*D'*C, -B*inv(Dg)*B';...
        C'*C+C'*D*inv(Dg)*D'*C, -A'-C'*D*inv(Dg)*B'];
    

w = linspace(0,20);
%preallocations
lam_max = 0; w_max =0;
for i = 1: length(w)
    %freq resp function and transpose
    P = C*inv(1j*w*eye(size(A))-A)*B+D;
    P_star = -B'*inv(1j*(2*pi*w(i))*eye(size(A))-A)*C'+D';
    if max(eig(P_star*P)) > lam_max 
        lam_max = max(eig(P_star*P));
        w_max = w(i);
    end
end
H_infn = lam_max;
end