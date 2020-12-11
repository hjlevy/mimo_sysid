function [Ass,Bss,Css,Dss] = model_generator(H,Htil,nmod)
%generating a model of size nmod based on larger hankel matrices

%taking the svd of the H100 matrix and decreasing to new state dim
[U,S,V] = svd(H);
Un = U(:,1:nmod); 
Sn = S(1:nmod,1:nmod); 
Vn = V(:,1:nmod);

%creating a new hankel matrix based on lower state dim
Hnew = Un*Sn*(Vn');

%Observability and Controllability
O = Un*diag(sqrt(diag(Sn)));
C = diag(sqrt(diag(Sn)))*Vn';

%C matrix in state model
Css = O(1:2,:);

%B matrix in state model
Bss = C(:,1:2);

%A matrix creation in state model
Oleft = inv(diag(sqrt(diag(Sn))))*Un';
Cright = Vn*inv(diag(sqrt(diag(Sn))));
Ass = Oleft*Htil*Cright;

Dss = zeros(2,2);
end