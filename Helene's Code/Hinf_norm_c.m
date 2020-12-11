function [H_infn,nu] = Hinf_norm_c(A,B,C,D,gam,tol)
%finding the Hinf norm and wmax of a continuous time system
I = eye(size(D'*D));
gamL = gam(1);
gamU = gam(2);

gamH = gamU;
while (gamH - gamL)/gamL > tol
    gamma = (gamL + gamH)/2;
    
    %building Aclp
    Dg = gamma^2*I-D'*D;
    Aclp = [A+B*inv(Dg)*D'*C, -B*inv(Dg)*B';...
            C'*C+C'*D*inv(Dg)*D'*C, -A'-C'*D*inv(Dg)*B'];
    lam = eig(Aclp);
    
    %if the minimum real part is close to zero (only imaginary part?)
    if min(abs(real(lam))) < 1.0e-5
        gamL = gamma;
    else
        gamH = gamma;
    end
end
H_infn = gamma;

%finding purely imaginary eigenvalue -> continuous time freq
nu = 0;
for k = 1:length(lam)
    %if almost no real part
    if  abs(real(lam(k))) < 1.0e-5
        nu = abs(imag(lam(k)));
    end
end

end