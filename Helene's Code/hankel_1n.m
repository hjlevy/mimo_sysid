function [Hnew,H] = hankel_1n(u,y,nr,nmod)
%finding Hankel Matrix for 1 input 1 output system
ind_off = find(u); %finds the index of the first non zero value

nc = nr;
H = zeros(nr,nc);
Htil = zeros(nr,nc);
for k = 1:nr
    %indexing based on offset when pulse hits
    ind = ind_off + k + [0:nc-1];
    
    H(k,:)=y(ind);
    Htil(k,:)=y(ind+1);
end

%taking the svd of the H100 matrix and decreasing to new state dim
[U,S,V] = svd(H);
Un = U(:,1:nmod); 
Sn = S(1:nmod,1:nmod); 
Vn = V(:,1:nmod);

%creating a new hankel matrix based on lower state dim
Hnew = Un*Sn*(Vn');
end
