clear
close
clc

A = 1;          % Amplitude
rb = 500;      % Bit rate
fs = 32*1000;   % Sampling frequency

Tb = 1/rb;      % Bit interval
Ts = 1/fs;      % Sampling interval


N  = 30;         % Number of bits to be transmitted


%% A distorted ( frequency selective) channe;

den = [1.00,-2.838,3.143,-1.709,0.458,-0.049];
num=0.1*[1,-1];
Ns = floor(Tb/Ts);

[H,w] = freqz(num,den,1024,fs); % Frequency response of the channel
figure(1)
plot(w/pi,20*log10(abs(H)),'LineWidth',2)
xlabel("Frequency")
ylabel("|H|")
title("Frequency response of the distorted channel")
set(gca,'fontsize',14)

%% Transmitted and Received signal

msg = randi([0,1],1,N);            % Binary data
[t,sig] = polar_nrz_lab5(msg,A,rb,fs);   % Polar Non return to zero scheme
distorted= filter(num,den,sig);     % Distorted signal
snr = 0;                            % Setting the SNR
received= awgn(distorted,snr);      % Received signal

figure(2)
plot(t,received,'LineWidth',1.5)
hold on;
plot(t,sig,'LineWidth',1.5)
hold on
xlabel("Time"), 
ylabel("Amplitude")
legend("Received","Transmitted")
title("Transmitted and received signal")
set(gca,'fontsize',14)

%% Matched filtering 

pulse = ones(1,Ns);                                     % Reference pulse
sig_filtered = filter(flip(pulse),1,received);          % Filtered signal
sig_filtered = sig_filtered/(A^2*Ns);                   % Normalizing

figure(3)
plot(t,sig_filtered,"LineWidth",2)
xlabel("Time"), 
ylabel("Amplitude")
title("Filtered")
set(gca,'fontsize',14)
%% Detection 

filtered = sig_filtered(Ns:Ns:end); % Sample at Tb
detected = filtered;
detected(filtered>=0)=1;            % Thresholding
detected(filtered<0)=0;
error = sum(msg~=detected)*(100/N);
disp(strcat("Error = ",num2str(error)," %"))
%% Eye diagram
eyediagram(sig_filtered,Ns);

