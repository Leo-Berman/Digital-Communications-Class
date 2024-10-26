clear all
close all
clc

%%
A = 1;              % Amplitude

rb = 1000;          % bit rate
Tb = 1/rb;          % bit interval     

fs = 100*rb;        % Sampling rate
Ts = 1/fs;          % Sampling Interval

M = 2;              % Numbe of symbols
N = 30;             % Number of bits to be transmitted

Ns = floor(Tb/Ts);  % Number of samples per bit

[t1,pnrz,t2,signal_prnz] = polar_nrz(A,rb,fs,N); % Polar Non return to zero scheme


figure(1)
stem(t1,pnrz);
xlabel('Time')
ylabel('Amplitude')

figure(2)
plot(t2,signal_prnz)
xlabel('Time')
ylabel('Amplitude')


%% Noisy Channel

snr_db = -10; % signal to noise power in db
received = awgn(signal_prnz,snr_db); % Passing the signal through a noisy channel
figure(3)
plot(t2,received)
xlabel('Time')
ylabel('Amplitude')

%% Matched filtering 
filtered_signal = zeros(1,N);

pulse = ones(1,Ns); % Reference pulse
% Performing matched filtering
for i = 1:length(signal_prnz)/Ns
    idx = (i-1)*Ns+1:(i-1)*Ns+Ns;
    filtered = conv(received(idx),pulse);
    filtered_signal(i) =  filtered(Ns);
end
filtered_signal = filtered_signal/(A^2*Ns);
%% Detection 
detected = filtered_signal;
detected(filtered_signal>=0)=1;
detected(filtered_signal<0)=-1;

figure(4)
stem(t1,detected);
xlabel('Time')
ylabel('Amplitude')

%% Performance 
error = sum(detected~=pnrz/A)/N



