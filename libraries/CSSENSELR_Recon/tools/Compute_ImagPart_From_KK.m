function  Imag_Data_f  = Compute_ImagPart_From_KK( Real_Data_f)
%COMPUTE_MINPHASE Summary of this function goes here

Imag_Data_f=-imag(hilbert(real((Real_Data_f)))); 


end

