function [t1,prz,t2,signal_upnrz] = SHANEM_unipolar_nrz(A,rb,fs,N)

    Tb = 1/rb; 
    Ts = 1/fs;
    M = 2;
     
    
    binary = randi([0,M-1],1,N);
    % prz = A*(2 * (binary - 0.5));
    prz=A.*binary;
    t1 = (0:N-1)*Tb; 
    
    upnrz = repelem(binary, 2);
    upnrz(2:2:end) = 0;
    
    signal_upnrz = repelem(upnrz, floor((Tb/2)/Ts));
    t2= (0:length(signal_upnrz)-1)*Ts;

end