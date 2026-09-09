function  Rephased_Data_trrr  = Compute_MinPhase( Data_trrr)
%COMPUTE_MINPHASE Summary of this function goes here
%MV=mean(fft(Data_trrr,[],1));
abs_Data_f=abs(fft(Data_trrr,[],1));

Rephased_Data_trrr=-ifft(abs_Data_f .* exp (-1j * imag (hilbert (log(abs_Data_f)))),[],1); 


end

