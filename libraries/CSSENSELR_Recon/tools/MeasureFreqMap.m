function [ FreqMap,CCoefMap,RefSpectrumShort ] = MeasureFreqMap(RefTimeSerie, Data_trrr,mrsiReconParams)

[~,MaxPPM_pt]=min(abs(mrsiReconParams.MaxPPM  - mrsiReconParams.ppm));
[~,MaxPPM_pt_B0Corr]=min(abs(mrsiReconParams.LRTGVModelParams.CorrB0Map_MaxPPM  - mrsiReconParams.ppm));
[~,MinPPM_pt_B0Corr]=min(abs(mrsiReconParams.LRTGVModelParams.CorrB0Map_MinPPM  - mrsiReconParams.ppm));

FreqRange=40;
FreqPrec=0.2;

Data_frrr=fft(Data_trrr,[],1);
Data_trrr=ifft(Data_frrr(MinPPM_pt_B0Corr:(end-MaxPPM_pt+MaxPPM_pt_B0Corr),:,:,:),[],1);

RefSpectrumShort=fft(RefTimeSerie,[],1);
RefSpectrumShort=RefSpectrumShort(MinPPM_pt_B0Corr:(end-MaxPPM_pt+MaxPPM_pt_B0Corr));
RefTimeSerie=ifft(RefSpectrumShort,[],1);

NbT=size(Data_trrr,1);
ImSize=size( mrsiReconParams.ImMask);
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;
Time=([0 :(NbT-1)]'/Fs+mrsiReconParams.AcqDelay);
PNSq_1rrr=(norm(RefTimeSerie)*sqrt(sum(abs(Data_trrr).^2,1) ) ).^2;

FreqArr=-FreqRange:FreqPrec:FreqRange;
CorCoef_rrrS= abs((RefTimeSerie' .*exp(2*pi*1i*Time*FreqArr).')* reshape(Data_trrr,NbT,[]) ).^2;
CorCoef_rrrS=reshape( transpose(CorCoef_rrrS),[ImSize numel(FreqArr)])./sqz(PNSq_1rrr);
FreqInd_111S(1,1,1,:)=FreqArr;

ThresCC_rrr=max(CorCoef_rrrS,[],4)/2;
pt_rrrS=(CorCoef_rrrS>ThresCC_rrr);
FreqMap=sum(CorCoef_rrrS.*pt_rrrS.*FreqInd_111S,4)./sum(CorCoef_rrrS.*pt_rrrS,4);
FreqMap(isnan(FreqMap))=0;
CCoefMap=max(CorCoef_rrrS,[],4);
CCoefMap(isnan(CCoefMap))=0;

FreqMap=FreqMap.*CCoefMap.*(CCoefMap>0.25);%/max(CCoefMap(mrsiReconParams.BrainMask>0));
FreqMap=(CCoefMap>0.25).*(FreqMap-mean(FreqMap( mrsiReconParams.ImMask>0 ) ) );% avoid freq. shifting

%smoothing of Maps
NX=size(Data_trrr,2);NY=size(Data_trrr,3);
sigma=mrsiReconParams.GaussianSigma;
[X,Y] = ndgrid(1:NX, 1:NY);
xc=floor(NX/2)+1;yc=floor(NY/2)+1;

exponent = -((X-xc).^2 + (Y-yc).^2)./(2*sigma^2);
Kernel = fftshift(fftshift(exp(exponent)/sum(exp(exponent(:))),1),2); 
Kernel=fft(fft(Kernel,[],1),[],2);

for s=1:size(FreqMap,3)
        FreqMap(:,:,s)=real(ifft(ifft(Kernel.*fft(fft(FreqMap(:,:,s),[],1),[],2),[],1),[],2));
end



FreqMap=FreqMap.*mrsiReconParams.ImMask;

end

