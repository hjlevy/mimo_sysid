function [H,Htil] = hankel_n(u,y1,y2,nr)
%finding Hankel Matrix for 2 input 2 output system
ind_off = find(u); %finds the index of the first non zero value

nc = nr;
H = zeros(2*nr,2*nc);
Htil = zeros(2*nr,2*nc);
for k = 1:nr
    %indexing based on offset when pulse hits
    ind = ind_off + k + [0:nc-1];
    
    %creating H based off response data with little h being 2x2 matrices
    %because 2 input 2 output system
    %  | h1 h2
    %  | h2 h3
    %(2k-1:2k) index needed in order to skip every 2 rows to make room for 
    %2x2 matrix
    H(2*k-1:2*k,1:2:end)=y1(:,ind);
    H(2*k-1:2*k,2:2:end)=y2(:,ind);
    
    %creating Htilda which is just H but shifted over 
    %  | h2 h3
    %  | h3 h4
    Htil(2*k-1:2*k,1:2:end)=y1(:,ind+1);
    Htil(2*k-1:2*k,2:2:end)=y2(:,ind+1);
end
end