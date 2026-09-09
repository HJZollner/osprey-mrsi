function  Rephased_Data_frrr  = Compute_MinPhase_f( Data_frrr)
%COMPUTE_MINPHASE Summary of this function goes here
%MV=mean(fft(Data_trrr,[],1));
abs_Data_f=abs(Data_frrr);

Rephased_Data_frrr=abs_Data_f .* exp (-1j * imag (hilbert (log(abs_Data_f)))); 


end

