%MAE 270A Project 
%Helene Levy
clear; close all; clc;

%% Loading and "massaging" initial data
load u1_impulse.mat
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi] = max(u1>0); %%% find index where pulse occurs
load u2_impulse.mat
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);

ts = 1/40; %%%% sample period
N = length(y11); %%%% length of data sets
t = [0:N-1]*ts - 1;

%% Plotting Original Data
% figure
% subplot(311)
% plot(t,u1,'b*','LineWidth',2)
% ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 1.1])
% 
% subplot(312)
% plot(t,y11,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 0.1])
% 
% 
% subplot(313)
% plot(t,y21,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% xlabel('second','FontSize',14)
% set(gca,'FontSize',14)
% axis([-0.2 2 -0.1 0.1])
% 
% figure
% subplot(311)
% plot(t,u2,'b*','LineWidth',2)
% ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 1.1])
% 
% subplot(312)
% plot(t,y12,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 0.1])
% 
% subplot(313)
% plot(t,y22,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% xlabel('second','FontSize',14)
% set(gca,'FontSize',14)
% axis([-0.2 2 -0.1 0.1])


%% Task 1 
%joining matrices of data representing same input
y1 = [y11; y21];
y2 = [y12; y22];

%plotting example hankel svds in prompt
ex = [20 40 80];
figure
for k = 1:length(ex)
    Hex = hankel_n(u1,y1,y2,ex(k));
    svd_vec = svd(Hex);
    k = 1:1:40;
    plot(k,svd_vec(k),'*');
    hold on;
    set(gca,'YScale', 'log');
    axis([0 40 10e-4 10e-1]);
    hold on;
end
xlabel('singular value index');
ylabel('Hankel singular value');
legend("H20","H40","H80");

%% Constructing Hankel Matrix and Plotting Singular Values
%constructing H100 matrix and plotting singular values
[H100,Htil] = hankel_n(u1,y1,y2,100);
svd_vec = svd(H100);

figure
k = 1:1:40;
plot(k,svd_vec(k),'*');
hold on;
set(gca,'YScale', 'log');
axis([0 40 10e-4 10e-1]);
xlabel('singular value index');
ylabel('Hankel singular value');
title('Singular Values of H100');

%% Simulating Impulse Response of Each State Model

%iterating through different state dims 
state_dim = [6,7,10,40];
An ={0}; Bn = {0}; Cn = {0}
for i = 1:length(state_dim)
    nmod = state_dim(i);
    %taking the svd of the H100 matrix and decreasing to new state dim
    [U,S,V] = svd(H100);
    U1 = U(:,1:nmod); 
    S1 = S(1:nmod,1:nmod); 
    V1 = V(:,1:nmod);

    %creating a new hankel matrix based on lower state dim
    Hnew = U1*S1*(V1');
    
    %Observability and Controllability
    O = U1*diag(sqrt(diag(S1)));
    C = diag(sqrt(diag(S1)))*V1';
    %C matrix in state model
    Css = O(1:2,:);
    Cn{i} = Css;
    %B matrix in state model
    Bss = C(:,1:2);
    Bn{i} = Bss;

    %A matrix creation in state model
    Oleft = inv(diag(sqrt(diag(S1))))*U1';
    Cright = V1*inv(diag(sqrt(diag(S1))));
    Ass = Oleft*Htil*Cright;
    An{i} = Ass;
    
    %checking stability 
    stab_check = max(abs(eig(Ass)));
    fprintf('Maximum eigenvalue of A for ns = %d, lambda = %4.3f\n',...
            nmod,stab_check);

    %simulating response to u1 = [1;0]
    h1 = zeros(2,100);
    x1 = Bss*[1;0]; %value x at k = 1
    for k = 1:100
        h1(:,k) = Css*x1;
        x1 = Ass*x1;
    end
    
    %simulating response to u2 = [0;1]
    h2 = zeros(2,100);
    x2 = Bss*[0;1]; %value x at k = 1
    for k = 1:100
        h2(:,k) = Css*x2;
        x2 = Ass*x2;
    end
    tsim = [1:100]*ts;

    %plotting each input/output response combo
    figure
    subplot(211)
    plot(tsim,h1(1,:),'*',t,y11,'*')
    grid on
    xlim([0 2]);
    ylabel('y1');
    title(sprintf('Impulse Response of u1 with ns = %4d',nmod));

    subplot(212)
    plot(tsim,h1(2,:),'*',t,y21,'*')
    grid on
    xlim([0 2]);
    ylabel('y2');
    xlabel('time (s)');
    
    figure
    subplot(211)
    plot(tsim,h2(1,:),'*',t,y12,'*')
    grid on
    xlim([0 2]);
    ylabel('y1');
    title(sprintf('Impulse Response pf u2 with ns = %4d',nmod));

    subplot(212)
    plot(tsim,h2(2,:),'*',t,y22,'*')
    grid on
    xlim([0 2]);
    ylabel('y2');
    xlabel('time (s)');
end
fprintf('State dimension 6 has an inferior reproduction of the data\n');

%% Creating and Plotting Model Frequency Response 

%plotting bode diagrams for different model order
N_w = 100;
w_span = linspace(0,20,N_w); %omega span in hertz
%iterating through different state dimensions
for k = 1:length(state_dim)
   
    %preallocating values to 0
    P11_mag = zeros(N_w,1);
    P12_mag = zeros(N_w,1);
    P21_mag = zeros(N_w,1);
    P22_mag = zeros(N_w,1);
    P11_ph = zeros(N_w,1);
    P12_ph = zeros(N_w,1);
    P21_ph = zeros(N_w,1);
    P22_ph = zeros(N_w,1);
    
    %function creating model based on each A,B,C matrix of state model
    %ex Cn = {{C6},{C7},{C10},{C40}}
    for i = 1:N_w
        %from equation of state C(e^jwtsI-A)^-1*B
        %making sure to convert w_span into rad/s for model analysis
        P = Cn{k}*inv(exp(1j*2*pi*w_span(i)*ts).*eye(size(An{k}))-An{k})...
            *Bn{k};
        
        %extracting freq model for y11
        P11_mag(i) = 20*log10(abs(P(1,1))); 
        P11_ph(i) = unwrap(angle(P(1,1)))*180/pi;
        
        %extracting freq model for y21
        P21_mag(i) = 20*log10(abs(P(2,1)));
        P21_ph(i) = unwrap(angle(P(2,1)))*180/pi;

        %extracting freq model for y12
        P12_mag(i) = 20*log10(abs(P(1,2)));
        P12_ph(i) = unwrap(angle(P(1,2)))*180/pi;
        
        %extracting freq model for y22
        P22_mag(i) = 20*log10(abs(P(2,2)));
        P22_ph(i) = unwrap(angle(P(2,2)))*180/pi;
    end

    %plotting model freq resp
    figure (11)
    subplot(211)
    plot(w_span,P11_mag);
    set(gca,'XScale', 'log');
    ylabel('Magnitude (db)');
    title('Bode Plot of System Response y1 to u1');
    grid on; hold on;

    subplot(212)
    plot(w_span,P11_ph);
    set(gca,'XScale', 'log');
    ylabel('Phase (deg)');
    xlabel('\omega (Hz)');
    grid on; hold on;

    figure (12)
    subplot(211)
    plot(w_span,P21_mag);
    set(gca,'XScale', 'log');
    ylabel('Magnitude (db)');
    title('Bode Plot of System Response y2 to u1');
    grid on; hold on;

    subplot(212)
    plot(w_span,P21_ph);
    set(gca,'XScale', 'log');
    ylabel('Phase (deg)');
    xlabel('\omega (Hz)');
    grid on; hold on;

    figure (13)
    subplot(211)
    plot(w_span,P12_mag);
    set(gca,'XScale', 'log');
    ylabel('Magnitude (db)');
    title('Bode Plot of System Response y1 to u2');
    grid on; hold on;

    subplot(212)
    plot(w_span,P12_ph);
    set(gca,'XScale', 'log');
    ylabel('Phase (deg)');
    xlabel('\omega (Hz)');
    grid on; hold on;

    figure (14)
    subplot(211)
    plot(w_span,P22_mag);
    set(gca,'XScale', 'log');
    ylabel('Magnitude (db)');
    title('Bode Plot of System Response y2 to u2');
    grid on; hold on;

    subplot(212)
    plot(w_span,P22_ph);
    set(gca,'XScale', 'log');
    ylabel('Phase (deg)');
    xlabel('\omega (Hz)');
    grid on; hold on;
end

%% Plotting Empirical Frequency Response 

figure (11);
subplot(211)
y11f = fft(y11)./fft(u1);
N = length(y11f);
om = [0:N-1]/(ts*N);
plot(om,20*log10(y11f));
set(gca,'XScale', 'log');
xlim([0 20]);
legend('ns = 6','ns = 7' ,'ns = 10' ,'ns = 40','empirical data',...
    'Location','southwest');
subplot(212)
plot(om,unwrap(angle((y11f)))*180/pi);
set(gca,'XScale', 'log');
xlim([0 20]);

figure (12);
subplot(211)
y21f = fft(y21)./fft(u1);
N = length(y21f);
om = [0:N-1]/(ts*N);
plot(om,20*log10(y21f));
set(gca,'XScale', 'log');
xlim([0 20]);
legend('ns = 6','ns = 7' ,'ns = 10' ,'ns = 40','empirical data',...
    'Location','southwest');
subplot(212)
plot(om,unwrap(angle((y21f)))*180/pi);
set(gca,'XScale', 'log');
xlim([0 20]);

figure (13);
subplot(211)
y12f = fft(y12)./fft(u2);
N = length(y12f);
om = [0:N-1]/(ts*N);
plot(om,20*log10(y12f));
set(gca,'XScale', 'log');
xlim([0 20]);
legend('ns = 6','ns = 7' ,'ns = 10' ,'ns = 40','empirical data',...
    'Location','southwest');
subplot(212)
plot(om,unwrap(angle((y12f)))*180/pi);
set(gca,'XScale', 'log');
xlim([0 20]);

figure (14);
subplot(211)
y22f = fft(y22)./fft(u2);
N = length(y22f);
om = [0:N-1]/(ts*N);
plot(om,20*log10(y22f));
set(gca,'XScale', 'log');
xlim([0 20]);
legend('ns = 6','ns = 7' ,'ns = 10' ,'ns = 40','empirical data',...
    'Location','southwest');
subplot(212)
plot(om,unwrap(angle((y22f)))*180/pi);
set(gca,'XScale', 'log');
xlim([0 20]);


