function [z,lam] = tzero_lam_gen(A,B,C,D,plot_bool)
%finding the transmission zeros and eigenvalues based on system matrices

ns = length(A); %state dimension
n = ns + size(B,2); %# columns of M matrix

%Constructing M matrix
%M = [A B; -C -D]
M = zeros(n,n);
M(1:ns,1:ns) = A;
M(1:ns,ns+1:end) = B;
M(ns+1:end,1:ns) = -C;
M(ns+1:end,ns+1:end) = -D;

%transmission zero calculation
z = eig(M,diag([ones(ns,1);zeros(n-ns,1)]));
z = sort(z,'descend'); %sorting from greatest magnitude to least 
inf_check = isinf(z);
z_ind = find(inf_check); %finding all the indices of "inf" values
z = z(z_ind(end)+1:end); %removing inf values

lam = eig(A);

if (plot_bool)
    figure;
    %plotting transmission zeros and eigenvalues
    plot(real(z),imag(z),'bo',real(lam),imag(lam),'rx'); 
    hold on;

    %plotting circle radius 1
    th = linspace(0,2*pi); r = 1;
    x_circ = r*cos(th);
    y_circ = r*sin(th);
    plot(x_circ,y_circ,'k--');
    xlabel('Re'); ylabel('Im');
    axis equal;grid on;
end
end