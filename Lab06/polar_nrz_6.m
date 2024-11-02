function [t2,signal_prz] = polar_nrz_6(binary, A,rb,fs)

Tb = 1/rb; 
Ts = 1/fs;

 


prz = A*(2 * (binary - 0.5));


signal_prz = repelem(prz, floor(Tb/Ts)); 
t2= (0:length(signal_prz)-1)*Ts;

end