%MAE 270A Project 
%Helene Levy
clear; close all; clc;

%plotting commands %change to false if don't want to plot
%task 1 plots
plot_resp = true; 
plot_fresp = false; 

%task 2 plots
tzero_plot = false;

%task 4 plots
plot_auto = false;
plot_cross = true;

%task 6 plots
hinf_plot = false;


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
% 

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
    j = 1:1:40;
    plot(j,svd_vec(j),'*');
    hold on;
    set(gca,'YScale', 'log');
    axis([0 40 10e-4 10e-1]);
    hold on;
end
xlabel('singular value index');
ylabel('Hankel singular value');
legend("H20","H40","H80");

%% Task 1 (cont) Constructing Hankel Matrix and Plotting Singular Values
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

    if(plot_resp)
        %plotting each input/output response combo
        figure
        subplot(211)
        plot(tsim,h1(1,:),'*',t,y11,'*');
        grid on
        xlim([0 2]);
        ylabel('y1');
        title(sprintf('Impulse Response of u1 with ns = %4d',nmod));
        legend('model','data','Location','southeast');

        subplot(212)
        plot(tsim,h1(2,:),'*',t,y21,'*');
        grid on
        xlim([0 2]);
        ylabel('y2');
        xlabel('time (s)');

        figure
        subplot(211)
        plot(tsim,h2(1,:),'*',t,y12,'*');
        grid on
        xlim([0 2]);
        ylabel('y1');
        title(sprintf('Impulse Response pf u2 with ns = %4d',nmod));
        legend('model','data','Location','southeast');

        subplot(212)
        plot(tsim,h2(2,:),'*',t,y22,'*');
        grid on
        xlim([0 2]);
        ylabel('y2');
        xlabel('time (s)');
    end
end
fprintf('State dimension 6 has an inferior reproduction of the data\n');

%% Task 1 (cont) Creating and Plotting Model Frequency Response 
%plotting bode diagrams for different model order
N_w = 200;
w_span = linspace(0,20,N_w); %omega span in hertz
%iterating through different state dimensions
for k = 1:length(state_dim)
    P = ss(An{k},Bn{k},Cn{k},Dn{k},ts);
    [P11_mag,P11_ph] = bode(P(1,1),2*pi*w_span);
    P11_mag = 20*log10(squeeze(P11_mag));
    P11_ph = squeeze(P11_ph);
    
    [P12_mag,P12_ph] = bode(P(1,2),2*pi*w_span);
    P12_mag = 20*log10(squeeze(P12_mag));
    P12_ph = squeeze(P12_ph);
    
    [P21_mag,P21_ph] = bode(P(2,1),2*pi*w_span);
    P21_mag = 20*log10(squeeze(P21_mag));
    P21_ph = squeeze(P21_ph);
    
    [P22_mag,P22_ph] = bode(P(2,2),2*pi*w_span);
    P22_mag = 20*log10(squeeze(P22_mag));
    P22_ph = squeeze(P22_ph);
    if plot_fresp
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
end


%% Task 1 (cont) Plotting Empirical Frequency Response 
if(plot_fresp)
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
end

%% Task 2: Transmission Zeros of MIMO Model
%2.1 finding transmission zeros of state model 7
ns = 7;
[A7,B7,C7,D7] = model_generator(H100,Htil,ns);
[z,lam] = tzero_lam_gen(A7,B7,C7,D7,tzero_plot);
title('Discrete Eigenvalues and Transmission Zeros');

%2.3 continuous time eigenvalues
lamc = log(lam)./ts;
%minimizing the omega by 2pik/ts for ones with imaginary part
%frequency bounded by [-pi/ts; pi/ts] rad/s which they already are

if(tzero_plot)
    figure;
    plot(real(lamc),imag(lamc),'rx'); 
    axis square ;grid on;
    title('Continuous Eigenvalues');
end

fprintf('From the set of lam_c we see there are two damped oscillators:\n')
disp('lam_c =');disp(lamc(1:4));
fprintf('Their natural frequences are %4.3f and %4.3f Hz \n',...
        imag(lamc(1))/(2*pi),imag(lamc(3))/(2*pi));
fprintf(['These frequencies look like the approximate cut-off ',... 
        'frequencies of the frequency responses\n\n']);
 
%% Task 2 (cont) Transmission Zeros for Each Channel 
%creating b1 b2 and c1 c2
b7_1 = B7(:,1); b7_2 = B7(:,2);
c7_1 = C7(1,:); c7_2 = C7(2,:);

%Channel (1,1) transmission zeros and eigenvalues
[z11,lam11] = tzero_lam_gen(A7,b7_1,c7_1,0,tzero_plot);
title(sprintf(['Channel (1,1): Discrete Eigenvalues and',...
    'Transmission Zeros for ns = %d'],ns));
fprintf(['3 state model: one oscillator and '...
        'one LP filter in channel (1,1)\n']);

%determining channel's singular values
H_11 = hankel_1n(u1,y11,100,3);
sig = svd(H_11);
disp('Singular Values for Channel (1,1):'); disp(sig(1:3));

%Channel (2,1) y2,u1
[z21,lam21] = tzero_lam_gen(A7,b7_1,c7_2,0,tzero_plot);
title(sprintf(['Channel (2,1): Discrete Eigenvalues and',...
    'Transmission Zeros for ns = %d'],ns));
fprintf('2 state model: one LP and one HP filter in channel (2,1)\n');

%determining channel's singular values
H_21 = hankel_1n(u1,y21,100,2);
sig = svd(H_21);
disp('Singular Values for Channel (2,1):'); disp(sig(1:2));

%Channel (1,2) y1,u2
[z12,lam12] = tzero_lam_gen(A7,b7_2,c7_1,0,tzero_plot);
title(sprintf(['Channel (1,2): Discrete Eigenvalues and',...
    'Transmission Zeros for ns = %d'],ns));
fprintf(['3 state model: one oscillator and one '...
        'LP filter in channel (1,2)\n']);
 
%determining channel's singular values
H_12 = hankel_1n(u2,y12,100,3);
sig = svd(H_12);
disp('Singular Values for Channel (1,2):'); disp(sig(1:3));

%Channel (2,2) y2,u2
[z22,lam22] = tzero_lam_gen(A7,b7_2,c7_2,0,tzero_plot);
title(sprintf(['Channel (2,2): Discrete Eigenvalues and',...
    'Transmission Zeros for ns = %d'],ns));
fprintf(['4 state model: one oscillator and one LP filter, '...
        'one HP filter in channel (2,2)\n']);

%determining channel's singular values
H_22 = hankel_1n(u2,y22,100,4);
sig = svd(H_22);
disp('Singular Values for Channel (2,2):'); disp(sig(1:4));

%% Task 2 (cont) Pole-Zero plot for nmod = 8
ns = 8;
[A8,B8,C8,D8] = model_generator(H100,Htil,ns);

%creating b1 b2 and c1 c2
b8_1 = B8(:,1); b8_2 = B8(:,2);
c8_1 = C8(1,:); c8_2 = C8(2,:);

%Channel (1,1) transmission zeros and eigenvalues
[z11,lam11] = tzero_lam_gen(A8,b8_1,c8_1,0,tzero_plot);
title(sprintf(['Channel (1,1): Discrete Eigenvalues and',...
    'Transmission Zeros for ns = %d'],ns));

%Channel (2,1) y2,u1
[z21,lam21] = tzero_lam_gen(A8,b8_1,c8_2,0,tzero_plot);
title(sprintf(['Channel (2,1): Discrete Eigenvalues and',...
    'Transmission Zeros for ns = %d'],ns));

%Channel (1,2) y1,u2
[z12,lam12] = tzero_lam_gen(A8,b8_2,c8_1,0,tzero_plot);
title(sprintf(['Channel (1,2): Discrete Eigenvalues and',...
    'Transmission Zeros for ns = %d'],ns));

%Channel (2,2) y2,u2
[z22,lam22] = tzero_lam_gen(A8,b8_2,c8_2,0,tzero_plot);
title(sprintf(['Channel (2,2): Discrete Eigenvalues and',...
    'Transmission Zeros for ns = %d'],ns));

%% Task 4: Impulse Response ID from White Noise Inputs
%parameters
load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1_rand = u_rand.Y(1).Data;
u2_rand = u_rand.Y(2).Data;
ts = 1/40;
N = length(y1);
t2 = [0:N-1]*ts - 1;

%remove any offsets from output data
y1 = y1 - mean(y1);
y2 = y2 - mean(y2);

u = [u1_rand;u2_rand];
y = [y1;y2];

%% Task 4 (cont) Autocorrelation of u and plots

%preallocations
% want k[-200,200] so give a 100 index offset on either end of data
q_min = 300; 
q_max = length(u1_rand)-300;
N = length(q_min:q_max)-1;

Ruu = zeros(2,401*2);
for i = 1:401
    %k in [-200, 200]
    k = i-201;
    
    %preallocation
    sum_Ruu = zeros(2,2);
    
    for q = q_min: q_max
        %summing each uk+q*uq^T over q set
        sum_Ruu = 1/(N)*(u(:,k+q)*u(:,q)')+sum_Ruu;
    end
    Ruu(:,2*i-1:2*i) = sum_Ruu;
end
Tau = linspace(-5,5,length(Ruu)/2);

%plotting autocorrelation graphs
if plot_auto
    figure;
    sgtitle('Autocorrelation of Inputs');
    subplot(221);
    %autocorr u1u1
    plot(Tau,Ruu(1,1:2:end)); hold on;
    axis([-5 5 -1 4]);
    xlabel('Lag Time, \tau');
    ylabel('Ruu');
    title('u1u1');

    %autocorr u1u2
    subplot(222);
    plot(Tau,Ruu(1,2:2:end)); hold on;
    axis([-5 5 -1 1]);
    xlabel('Lag Time, \tau');
    ylabel('Ruu');
    title('u1u2');

    %autocorr u2u1
    subplot(223);
    plot(Tau,Ruu(2,1:2:end)); hold on;
    axis([-5 5 -1 1]);
    xlabel('Lag Time, \tau');
    ylabel('Ruu');
    title('u2u1');

    %autocorr u2u2
    subplot(224);
    plot(Tau,Ruu(2,2:2:end)); hold on;
    axis([-5 5 -1 4]);
    xlabel('Lag Time, \tau');
    ylabel('Ruu');
    title('u2u2');
end
%variance of input matrix (Ruu[0]);
V = Ruu(:,401:402);
disp('Ruu[0]=');disp(V);
fprintf('This indicates a variance of ~4 for each input channel\n');

%% Task 4 (cont) Cross correlation of u,y and plots

%preallocations
% want k[-200,200] so give a 100 index offset on either end of data
q_min = 300; 
q_max = length(u1_rand)-300;
N = length(q_min:q_max)-1;

Ruu = zeros(2,401*2);
%lag time corresponds to T = k*ts -> k = Tbound/ts
kmin = -0.2/ts; kmax= 2/ts;
for i = 1:kmax-kmin+1
    %k in [-200, 200]
    %k is [-0.2/ts,2/ts]
    k = i+(kmin-1);
    
    %preallocation
    sum_Ryu = zeros(2,2);
    
    for q = q_min: q_max
        %summing each uk+q*uq^T over q set
        sum_Ryu = 1/(N)*(y(:,k+q)*u(:,q)')+sum_Ryu;
    end
    Ryu(:,2*i-1:2*i) = sum_Ryu;
end
%normalizing columns by variance of u1 data
Ryu(:,1:2:end) = Ryu(:,1:2:end)/(V(1,1));

%normalizing columns by variance of u2 data
Ryu(:,2:2:end) = Ryu(:,2:2:end)/(V(2,2));
Tau = linspace(-0.2,2,length(Ryu)/2);

%plotting cross correlation graphs
if plot_cross
    figure;
    subplot(211);
    %crosscorr y1u1
    plot(Tau,Ryu(1,1:2:end),'b*',t,y11,'r*'); 
    axis([0 1 -0.1 0.1]);
    xlabel('time');
    ylabel('y11');
    legend('Ryu','pulse exp');
    title('Cross Correlation y1,u1');
    
    %crosscorr y2u1
    subplot(212);
    plot(Tau,Ryu(2,1:2:end),'b*',t,y21,'r*'); 
    axis([0 1 -0.1 0.1]);
    xlabel('time');
    ylabel('y21');
    legend('Ryu','pulse exp');
    title('Cross Correlation y2,u1');


    %crosscorr y1u2
    figure;
    subplot(211);
    plot(Tau,Ryu(1,2:2:end),'b*',t,y12,'r*'); 
    axis([0 1 -0.1 0.1]);
    xlabel('time');
    ylabel('y12');
    legend('Ryu','pulse exp');
    title('Cross Correlation y1,u2');

    %crosscorr y2u2
    subplot(212);
    plot(Tau,Ryu(2,2:2:end),'b*',t,y22,'r*');
    axis([0 1 -0.1 0.1]);
    xlabel('time');
    ylabel('y22');
    legend('Ryu','pulse exp');
    title('Cross Correlation y2,u2');
end

%% Task 4 (cont) Autocorrelation of y (scaled)
y_scale = y/2;
%preallocations
% want k[-200,200] so give a 100 index offset on either end of data
q_min = 300; 
q_max = length(y1)-300;
N = length(q_min:q_max)-1;

Ryy = zeros(2,401*2);
for i = 1:401
    %k in [-200, 200]
    k = i-201;
    
    %preallocation
    sum_Ryy = zeros(2,2);
    
    for q = q_min: q_max
        %summing each uk+q*uq^T over q set
        sum_Ryy = 1/(N)*(y_scale(:,k+q)*y_scale(:,q)')+sum_Ryy;
    end
    Ryy(:,2*i-1:2*i) = sum_Ryy;
end
Ryy_0 = Ryy(:,401:402);

%% Task 5: H2 norm
%RMS value of output when u is zero mean unit variance noise
y_rms = sqrt(trace(Ryy_0));
fprintf('RMS value of scaled output y = %4.4f\n',y_rms);

%||P||H2 from 7 state model
ns = 7;
[A7,B7,C7,D7] = model_generator(H100,Htil,ns);

P_norm9 = 0;
for k = 1: length(u1)
    P_norm9 = P_norm9 + trace(B7'*(A7')^k*(C7')*C7*(A7^k)*B7);
end
P_norm9 = sqrt(P_norm9);
fprintf('||P||H2 using State Space Model Eqn (9) = %4.4f\n',P_norm9);

P_norm10 = 0;
for k = 1: length(u1)
    P_norm10 = P_norm10 + trace(C7*(A7^k)*B7*(B7')*(A7')^k*C7');
end
P_norm10 = sqrt(P_norm10);
fprintf('||P||H2 using State Space Model Eqn (10) = %4.4f\n',P_norm10);

%||P||H2 from discrete time pulse response
y_imp = [y11 y12; y21 y22];
P_norm8 = 0; 
for k = 1: length(u1)
    P_norm8 = P_norm8 + trace(y_imp(:,2*k-1:2*k)'*y_imp(:,2*k-1:2*k));
end
P_norm8 = sqrt(P_norm8);
fprintf('||P||H2 using discrete time pulse response y Eqn (8) = %4.4f\n',...
        P_norm8);

%% Task 6: H infinity norm

plot_6fresp = true;
%H infinity norm calculations
svd_vec = svd(H100); ns = 7;
%recommended bounds for gamma found online
gam = [svd_vec(1),2*sum(svd_vec)];
tol = 0.01;
[Hinf, wmax] = Hinf_norm_d(A7,B7,C7,D7,gam,tol,ts);
fprintf('Hinf norm = %4.3f = %4.3f db\n',Hinf,20*log10(Hinf));
fprintf('where max singular value achieves largest value, w = %4.3f Hz\n', ...
    wmax/(2*pi));

%frequency response plot using model
P = ss(A7,B7,C7,D7,ts);
[SV,W] = sigma(P,2*pi*w_span); hold on;
if hinf_plot
    figure;
    plot(W/(2*pi),20*log10(SV(1,:)));  hold on;
    plot(W/(2*pi),20*log10(SV(2,:))); hold on;
end

%frequency response plot using pulse response data
y11f = fft(y11)./fft(u1);
y21f = fft(y21)./fft(u1);
y12f = fft(y11)./fft(u2);
y22f = fft(y22)./fft(u2);

%forming frequency model values for each corresponding frequency
N = length(y11f);
om = [0:N-1]/(ts*N);
sig = zeros(2,length(om));
for i = 1:length(om)
    H = [y11f(i) y12f(i); y21f(i) y22f(i)];
    %taking singular values of each H
    sig(:,i) = svd(H);
end

%plotting frequency vs. singular values of pulse response data
if hinf_plot
    plot(om,20*log10(sig(1,:)));  hold on;
    plot(om,20*log10(sig(2,:))); hold on;
    plot(wmax/(2*pi),20*log10(Hinf),'k*'); hold on;
    legend('y1 model','y2 model','y1 data','y2 data','Hinf','Location','southwest');
    axis([0,20,-55, -5])
    set(gca,'XScale', 'log');
    xlabel('\omega (Hz)');
    ylabel('Singular Values (dB)');
    title('Singular Values of Frequency Response');
end
