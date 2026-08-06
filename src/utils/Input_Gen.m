%% Input Signal generator code - JP 29/-9/23

%clc; close all; clear all;
L = 5000; %number of time step in CFD
U_inf=0.01460734693878; %Inlet velocity
mu = 0;
nu=0.05;
fshed=0.0032; %Hz
fc = 0.75*fshed;
Ts = 0.4; %time step
t0 = 10400; %Initial time
t = t0:Ts:Ts*(10000-1)+t0;
fs =1/Ts;
U=normrnd(mu,nu,L,1); %Generate random signal
O=3; %Order of butter
[b,a]=butter(O,fc/(fs/2)); %filter
U=[zeros(5000,1); U];
Cmu = (2*(U.^2*0.05)/(U_inf^2));

figure(1)
freqz(b,a,[],fs) %plot filter char
subplot(2,1,1)
%xlim([1 40]) 
%ylim([-100 20])
x=filter(b,a,U);


figure(2)
L = size(x,1);             % Length of signal
Y = fft(x);
plot(fs/L*(0:L-1),abs(Y),"LineWidth",3)
%pwelch(x, hamming(size(x,1)), [ ], [ ], Fs)
%xlim([0 0.08])
title("Complex Magnitude of fft Spectrum")
xlabel("f (Hz)")
ylabel("|fft(X)|")
U=x;


figure(4)
plot(t,Cmu(:,1), LineWidth=2); %xlim([0 10e-4]) %plot input signal
xlabel('Time (s)');
ylabel('Momentum Coefficient')
grid on
p = 0.5.*U.^2.*1.225;
rms(x) %check rms

%p=[zeros(2000,1); p];
U=[t' U];

figure(3)
plot(U(:,1),U(:,2)/13.27)
%P=[T p];
writematrix(U,'U.txt','Delimiter','space') %output signal
%writematrix(P,'P.txt','Delimiter','space') %output signal
%%
% theta=85:0.1:95;%85:0.1:95;
% theta=theta*pi/180;
% x=0.05*cos(theta);
% w=10;
% L=0.32;
% rho=1.225;
% mdes=1400;
% vy=((9*((mdes*rho)/60000))/(rho*L*0.05)).*cos(18*theta-9*pi); %3/320 LPM to velocity profile
% figure
% plot(x,vy, LineWidth=2)
% xlabel('x-coord');
% ylabel('Velocity, m/s')
% grid on