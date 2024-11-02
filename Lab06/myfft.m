function [f,psd] = myfft(signal,N,fs)

X = (1/N)*fftshift(fft(signal,N));
df = fs/N;
sampleIndex = -N/2:1:N/2-1;
f = sampleIndex*df;

% magnitude spectrum
mag = abs(X);

% Phase spectrum
phase = atan2(imag(X),real(X));
X2 = X;
tol = max(abs(X))/10000;
X2(abs(X)<tol) = 0;
phase = rad2deg(atan2(imag(X2),real(X2)));

% PSD
psd = mag.^2;

end