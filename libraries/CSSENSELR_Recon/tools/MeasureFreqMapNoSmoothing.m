function [ FreqShift_rrr,MaxCCoef_rrr,RefSpectrumShort ] = MeasureFreqMapNoSmoothing(RefTimeSerie, Data_trrr,Mask,mrsiReconParams)

NbT=size(Data_trrr,1);
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;

FreqRange=mrsiReconParams.WaterRefFreqSearchRange*2; %hz
FreqPrec=0.5; %hz

%[~,MaxPPM_pt]=min(abs(mrsiReconParams.MaxPPM - mrsiReconParams.ppm));
%[~,MaxPPM_pt_B0Corr]=min(abs(mrsiReconParams.LRTGVModelParams.CorrB0Map_MaxPPM  - mrsiReconParams.ppm));
%[~,MinPPM_pt_B0Corr]=min(abs(mrsiReconParams.LRTGVModelParams.CorrB0Map_MinPPM  - mrsiReconParams.ppm));

%if MaxPPM_pt<MaxPPM_pt_B0Corr
%	MaxPPM_pt_B0Corr = MaxPPM_pt;
%end

Cutoff=round((FreqRange/2)*NbT/Fs);


Data_frrr=fft(Data_trrr,[],1);
%Data_trrr=ifft(Data_frrr(MinPPM_pt_B0Corr:(end-MaxPPM_pt+MaxPPM_pt_B0Corr),:,:,:),[],1);
Data_trrr=ifft(Data_frrr(Cutoff:(end-Cutoff),:,:,:),[],1);

RefSpectrumShort=fft(RefTimeSerie,[],1);
%RefSpectrumShort=RefSpectrumShort(MinPPM_pt_B0Corr:(end-MaxPPM_pt+MaxPPM_pt_B0Corr));
RefSpectrumShort=RefSpectrumShort(Cutoff:(end-Cutoff));
RefTimeSerie=ifft(RefSpectrumShort,[],1);

 %Recompute dimension and samplerate after cutoff
NbT=size(Data_trrr,1);
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;
Time=([0 :(NbT-1)]'/(Fs));

PNSq_1rrr=(norm(RefTimeSerie)*sqrt(sum(abs(Data_trrr).^2,1) ) ).^2;
shfData_trrr=0*Data_trrr;
FreqArr=-FreqRange:FreqPrec:FreqRange;

CorCoef_rrrS=zeros([size( mrsiReconParams.BrainMask) numel(FreqArr)],class(Data_trrr));
FreqInd_111S=zeros([1 1 1 numel(FreqArr)],class(Data_trrr));
ind=0;
for FreqS=FreqArr
    ind=ind+1;
    shfData_trrr=(Data_trrr.*exp(2*pi*1i*Time*FreqS));
    CorCoef_rrrS(:,:,:,ind)=reshape( abs(RefTimeSerie' * reshape(shfData_trrr,NbT,[]) ).^2,size( PNSq_1rrr) )./PNSq_1rrr;
    FreqInd_111S(1,1,1,ind)=FreqS;
end

ThresCC_rrr=max(CorCoef_rrrS,[],4)/2;
pt_rrrS=(CorCoef_rrrS-ThresCC_rrr)>0;
FreqShift_rrr=sum(CorCoef_rrrS.*pt_rrrS.*FreqInd_111S,4)./sum(CorCoef_rrrS.*pt_rrrS,4);
MaxCCoef_rrr=max(CorCoef_rrrS,[],4);

clear pt_rrrS CorCoef_rrrS FreqInd_111S ThresCC_rrr

MaxCCoef_rrr(isnan(MaxCCoef_rrr(:)))=0;
FreqShift_rrr(isnan(FreqShift_rrr(:)))=0;

MaxCCoef_rrr = Mask.*MaxCCoef_rrr;
%Weighting
FreqShift_rrr=Mask.*FreqShift_rrr.*MaxCCoef_rrr.*(MaxCCoef_rrr>0.25);




end

