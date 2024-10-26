clear all
close all
clc

%% Data stream generation

% Parameters
A = 1; 
rb = 1000;      % Bit rate
Tb = 1/rb;      % Bit interval
fs = 100*rb;    % Sampling freq
Ts = 1/fs;      % Sampling interval

M = 2;                  % Number of symbols
Nb = ceil(log2(M));     % Number of bits per symbol

N = 32;
t1 = (0:N-1)*Tb;
binary = A*randi([0,M-1],1,N);

figure(1)
stem(t1,binary)
xlabel('Time')
ylabel('Binary Data')
title('Binary Data')

%% Increasing sampling rate (Also Unipolar nonreturn to zero)
signal = repelem(binary,floor(Tb/Ts));
t = (0:length(signal)-1)*Ts;
figure(2)
plot(t,signal)
xlabel('Time')
ylabel('UPNRZ')
title('UPNRZ')

[f,psd]=myfft(signal,length(signal),fs);
psd_dbw = 10*log10(psd);

figure(3)
plot(f(f>=0),psd_dbw(f>=0));
xlabel('dBW')
ylabel('Frequency')
title('Power Spectral Density of UPNRZ ')



scope1 = dsp.SpectrumAnalyzer();
scope1.SampleRate = fs;
scope1.PlotAsTwoSidedSpectrum = false;
scope1.SpectrumUnits = "dBW";
scope1(signal');
scope1.Name = 'Unipolar NRZ';
release(scope1);

%% Polar NRZ
bnrz = 1*(2 * (binary - 0.5));signal_bnrz = repelem(bnrz, floor(Tb/Ts)); 

figure(4)
plot(t,signal_bnrz)
xlim([0 .02])
xlabel('Time')
ylabel('PNRZ')
title('PNRZ')
% bb = signal;
% bb(signal==0)=-1;
% figure(5)
% plot(t,bb)

[f,psd]=myfft(signal_bnrz,length(signal),fs);
psd_dbw = 10*log10(psd);


figure(5)
plot(f(f>=0),psd_dbw(f>=0));
xlabel('dBW')
ylabel('Frequency')
title('Power Spectral Density of PNRZ ')

scope2 = dsp.SpectrumAnalyzer();
scope2.SampleRate = fs;
scope2.PlotAsTwoSidedSpectrum = false;
scope2.SpectrumUnits = "dBW";
scope2(signal_bnrz');
scope2.Name = 'Polar RZ';
release(scope2);

%% Unipolar return to zero
urz = repelem(binary, 2);
urz(2:2:end) = 0;
signal_urz = repelem(urz, floor((Tb/2)/Ts)); 
figure(6)
plot(t,signal_urz)
xlim([0 .02])
xlabel('Time')
ylabel('UPRZ')
title('UPRZ')

[f,psd]=myfft(signal_urz,length(signal),fs);
psd_dbw = 10*log10(psd);

figure(7)
plot(f(f>=0),psd_dbw(f>=0));
xlabel('dBW')
ylabel('Frequency')
title('Power Spectral Density of UPRZ ')

scope3 = dsp.SpectrumAnalyzer();
scope3.SampleRate = fs;
scope3.PlotAsTwoSidedSpectrum = false;
scope3.SpectrumUnits = "dBW";
scope3(signal_urz');
scope3.Name = 'Unipolar RZ';
release(scope3);

%% Manchester
% Manchester
manchester = repelem(binary, 2);
manchester_zeros = 0 * manchester(manchester == 0); 
manchester_zeros(1:2:end) = -1; 
manchester_zeros(2:2:end) = 1;

manchester_ones = 0 * manchester(manchester == 1); 
manchester_ones(1:2:end) = 1; 
manchester_ones(2:2:end) = -1;

man = manchester;
man(manchester == 0) = manchester_zeros;
man(manchester == 1) = manchester_ones;
signal_man = repelem(man, floor((Tb/2)/Ts));

figure(8)
plot(t,signal_man)
xlim([0 .02])
xlabel('Time')
ylabel('Manchester Coding')
title('Manchester Coding')

[f,psd]=myfft(signal_man,length(signal),fs);
psd_dbw = 10*log10(psd);

figure(9)
plot(f(f>=0),psd_dbw(f>=0));
xlabel('dBW')
ylabel('Frequency')
title('Power Spectral Density of manchester coding ')

scope4 = dsp.SpectrumAnalyzer();
scope4.SampleRate = fs;
scope4.PlotAsTwoSidedSpectrum = false;
scope4.SpectrumUnits = "dBW";
scope4(signal_man');
scope4.Name = 'Manchester';
release(scope4);
%% Differential 

dunrz = zeros(size(binary));

ref = 1;
for c =1:length(binary)
    dunrz(c) = not(xor(ref,binary(c)));
    ref = dunrz(c);


end

signal_dunrz = repelem(dunrz, floor(Tb/Ts));       % Change indices to the new sampling frequency

figure(10)
plot(t, signal_dunrz); 
xlim([0 .02])
xlabel('Time')
ylabel('Differential')
title('Differential Encoding')

[f,psd]=myfft(signal_dunrz,length(signal),fs);
psd_dbw = 10*log10(psd);

figure(11)
plot(f(f>=0),psd_dbw(f>=0));
xlabel('dBW')
ylabel('Frequency')
title('Power Spectral Density of Differential Encoding ')

scope5 = dsp.SpectrumAnalyzer();
scope5.SampleRate = fs;
scope5.PlotAsTwoSidedSpectrum = false;
scope5.SpectrumUnits = "dBW";
scope5(signal_dunrz');
scope5.Name = 'Differential';
release(scope5);


