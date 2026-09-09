function Filtered_Data_trrr=Rephase_Data(Data_trrr,reconResults)
Filtered_Data_trrr=zeros(size(Data_trrr));
Spectrum=zeros(size(Data_trrr,1),1);
for a=1:size(Data_trrr,2)
    for b=1:size(Data_trrr,3)
        for c=1:size(Data_trrr,4)
            if(reconResults.BrainMask(a,b,c))
                %  Filtered_Data_trr(:,a,b)=KKRecurs(Data_trrr(:,a,b));
                Spectrum=fft(Compute_MinPhase(Data_trrr(:,a,b,c)));
                Spectrum(reconResults.MinPPM_pt:reconResults.MaxPPM_pt)=Rephase_poly(Spectrum(reconResults.MinPPM_pt:reconResults.MaxPPM_pt));
                Filtered_Data_trrr(:,a,b,c)=ifft(Spectrum);
            end
        end
    end
end
end
