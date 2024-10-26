function [t,signal_man] = manchester_sig(binary,A,rb,fs)

Tb = 1/rb; 
Ts = 1/fs;

 


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
prz1 = A*(2 * (binary - 0.5));


t = (0:length(signal_man)-1)*Ts;

end