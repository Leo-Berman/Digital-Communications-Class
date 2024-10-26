function [t1,prz,t2,signal_prz] = polar_nrz(A,rb,fs,N)

Tb = 1/rb; 
Ts = 1/fs;
M = 2;
 

binary = randi([0,M-1],1,N);
prz = A*(2 * (binary - 0.5));
t1 = (0:N-1)*Tb; 

signal_prz = repelem(prz, floor(Tb/Ts)); 
t2= (0:length(signal_prz)-1)*Ts;

end