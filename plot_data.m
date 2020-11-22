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
An ={0}; Bn = {0}; Cn = {0}; Dn = {0};
for i = 1:length(state_dim)
    nmod = state_dim(i);
    
    [An{i},Bn{i},Cn{i},Dn{i}] = model_generator(H100,Htil,nmod);

    %simulating response to u1 = [1;0]
    h1 = zeros(2,100);
    x1 = Bn{i}*[1;0]; %value x at k = 1
    for k = 1:100
        h1(:,k) = Cn{i}*x1;
        x1 = An{i}*x1;
    end
    
    %simulating response to u2 = [0;1]
    h2 = zeros(2,100);
    x2 = Bn{i}*[0;1]; %value x at k = 1
    for k = 1:100
        h2(:,k) = Cn{i}*x2;
        x2 = An{i}*x2;
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

%% Task 2
%2.1 finding transmission zeros of state model 7
ns = 7;
A7 = An{2}; C7 = Cn{2}; B7 = Bn{2}; D7 = Dn{2};

%Constructing M matrix
%M = [A7 B7; -C7 -D7]
M = zeros(9,9);
M(1:ns,1:ns) = A7;
M(1:ns,ns+1:end) = B7;
M(ns+1:end,1:ns) = -C7;
M(ns+1:end,ns+1:end) = -D7;

%finding transmission zeros from M
%Mv = z[I 0;0 0]v
z = eig(M,diag([ones(7,1);zeros(2,1)]));
z = sort(z,'descend'); %sorting from greatest magnitude to least 
z = z(5:9,1); %removing infinity values

%2.2 plotting eigenvalues and tzeros
%eigenvalues of the system
lam = eig(A7);
figure;
plot(real(z),imag(z),'bo',real(lam),imag(lam),'rx'); hold on;
title('Discrete Eigenvalues and Transmission Zeros');

%plotting circle radius 1
th = linspace(0,2*pi); r = 1;
x_circ = r*cos(th);
y_circ = r*sin(th);
plot(x_circ,y_circ,'k--');
axis equal;grid on;

%2.3 continuous time eigenvalues
lamc = log(lam)./ts;
%minimizing the omega by 2pik/ts for ones with imaginary part
%frequency bounded by [-pi/ts; pi/ts] rad/s which they already are

figure;
plot(real(lamc),imag(lamc),'rx'); hold on;
axis square ;grid on;
title('Continuous Eigenvalues');

fprintf('From the set of lam_c we see there are two damped oscillators:\n')
disp('lam_c =');disp(lamc(1:4));
fprintf('Their natural frequences are %4.3f and %4.3f Hz \n',...
        imag(lamc(1))/(2*pi),imag(lamc(3))/(2*pi));
fprintf(['These frequencies look like the approximate cut-off',... 
        'frequencies of the frequency responses\n']);
 
%% Transmission Zeros for Each Channel (2.4)
%creating b1 b2 and c1 c2
b7_1 = B7(:,1); b7_2 = B7(:,2);
c7_1 = C7(1,:); c7_2 = C7(2,:);

%Channel 1
M7_1 = zeros(8,8);
M7_1(1:ns,1:ns) = A7;
M7_1(1:ns,ns+1:end) = b7_1;
M7_1(ns+1:end,1:ns) = -c7_1;
M7_1(ns+1:end,ns+1:end) = 0;

%trasmission zeros of channel 1
z = eig(M7_1,diag([ones(7,1);0]));
z = sort(z,'descend'); %sorting from greatest magnitude to least 
z = z(3:end); %removing inf values

figure;
plot(real(z),imag(z),'bo',real(lam),imag(lam),'rx'); hold on;
title('Channel 1: Discrete Eigenvalues and Transmission Zeros nmod=7');

%analyzing pole, zero cancellations from figure
%oscillators have real part and imaginary part
%remember e^(sig+jw)t = e^(sigt)(sinwt+coswt)
fprintf('Channel 1 has one oscillator and one LP Filter\n')
fprintf('This means this channel can be represented by 2 states.\n')


%Channel 2
M7_2 = zeros(8,8);
M7_2(1:ns,1:ns) = A7;
M7_2(1:ns,ns+1:end) = b7_2;
M7_2(ns+1:end,1:ns) = -c7_2;
M7_2(ns+1:end,ns+1:end) = 0;

%transmission zeros of channel 2
z = eig(M7_2,diag([ones(7,1);0]));
z = sort(z,'descend'); %sorting from greatest magnitude to least 
z = z(3:end); %removing inf values

figure;
plot(real(z),imag(z),'bo',real(lam),imag(lam),'rx'); hold on;
title('Channel 2: Discrete Eigenvalues and Transmission Zeros nmod=7');

%analyzing pole, zero cancellations from figure
fprintf('Channel 2 has one oscillator and two LP Filter.\n')
fprintf('This means this channel can be represented by 3 states.\n')

%#4 WHAT ARE HANKEL SINGULAR VALUES OF EACH CHANNEL ?

%% Pole-Zero plot for nmod = 8
ns = 8;
[A8,B8,C8,D8] = model_generator(H100,Htil,ns);

%channel 1 M
M8 = zeros(9,9);
M8(1:ns,1:ns) = A8;
M8(1:ns,ns+1:end) = B8(:,1);
M8(ns+1:end,1:ns) = -C8(1,:);
M8(ns+1:end,ns+1:end) = 0;

%transmission zeros of channel 1
z = eig(M8,diag([ones(8,1);0]));
z = sort(z,'descend'); %sorting from greatest magnitude to least 
z = z(3:end); %removing inf values

lam = eig(A8);
figure;
plot(real(z),imag(z),'bo',real(lam),imag(lam),'rx'); hold on;
title('Channel 1: Discrete Eigenvalues and Transmission Zeros nmod=8');
%pole zero cancellation analysis
fprintf('Channel 1 of model dimension 8 only requires 2 states');

%channel 2 M
M8 = zeros(9,9);
M8(1:ns,1:ns) = A8;
M8(1:ns,ns+1:end) = B8(:,2);
M8(ns+1:end,1:ns) = -C8(2,:);
M8(ns+1:end,ns+1:end) = 0;

%transmission zeros of channel 2
z = eig(M8,diag([ones(8,1);0]));
z = sort(z,'descend'); %sorting from greatest magnitude to least 
z = z(3:end); %removing inf values

figure;
plot(real(z),imag(z),'bo',real(lam),imag(lam),'rx'); hold on;
title('Channel 2: Discrete Eigenvalues and Transmission Zeros nmod=8');
%pole zero cancellation analysis
fprintf('Channel 2 of model dimension 8 only requires 3 states');