function [H,Htil] = hankel_1n(u,y,nr)
%finding Hankel Matrix for 2 input 2 output system
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
end