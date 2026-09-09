function [MRSIData_ctkkk ]=WaterSuppressionProj(MRSIData_ctkkk,mrProt,FiltParam,kmask)
NSample=1000;
[NbCoil,NbT,Nr,Nc,Ns]=size(MRSIData_ctkkk);
temp_ckkk_t = reshape(permute(MRSIData_ctkkk,[1,3,4,5,2]),[],NbT);

clear MRSIData_ctkkk

Energy_ckkk=sum(abs(temp_ckkk_t).^2,2);

[~, I] =sort(Energy_ckkk,'descend');
WaterFree_st=zeros(NSample,NbT);
for samp=1:NSample
    SampleSpectra_st(samp,:) =  Fast_HSVD_Filter(transpose(temp_ckkk_t(I(samp),:)),1.0/mrProt.DwellTime,32,FiltParam.Water_minFreq,FiltParam.Water_maxFreq);
end

Water_st = temp_ckkk_t(I(1:NSample),:) - SampleSpectra_st;
[~,Sorig,Vorig] = svd(Water_st,0);

WaterProj=Vorig(:,1:FiltParam.Comp)*Vorig(:,1:FiltParam.Comp)';

IOP=diag(ones(NbT,1));

MRSIData_ctkkk =	temp_ckkk_t*(IOP-WaterProj);

MRSIData_ctkkk=permute(reshape(MRSIData_ctkkk,[NbCoil,Nr,Nc,Ns,NbT]),[1,5,2,3,4]);

% figure;
% TEST =  Fast_HSVD_Filter(permute(temp_ckkk_t(1,:),[2,1]),1.0/mrProt.DwellTime,128,FiltParam.Water_minFreq,FiltParam.Water_maxFreq);
% hold on 
% plot(abs(fft(temp_ckkk_t(1,:),[],2)))
% plot(abs(fft(MRSIData_ctkkk(1,:,1,1,1),[],2)))
% plot(abs(fft(TEST)))


end
