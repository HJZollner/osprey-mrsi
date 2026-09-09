function [MRSIData_ctkkk ]=WaterSuppressionProjHSVD(MRSIData_ctkkk,mrProt,FiltParam,kmask)

[NbCoil,NbT,Nr,Nc,Ns]=size(MRSIData_ctkkk);
temp_ckkk_t = reshape(permute(MRSIData_ctkkk,[1,3,4,5,2]),[],NbT);

clear MRSIData_ctkkk

[~,Sorig,Vorig] = svd(temp_ckkk_t,0);
NB_Comp_FullWater=floor(NbT/4); % constrained by the Hankel Matrix in HSVD
if NB_Comp_FullWater>size(Vorig,2)
   NB_Comp_FullWater=size(Vorig,2);
end
V_tc=conj(Vorig(:,1:NB_Comp_FullWater)*Sorig(1:NB_Comp_FullWater,1:NB_Comp_FullWater));
FV_tc = 0*V_tc;
for c=1:size(V_tc,2);
    
    FV_tc(:,c) =  Fast_HSVD_Filter(V_tc(:,c),1.0/mrProt.DwellTime,FiltParam.Comp,FiltParam.Water_minFreq,FiltParam.Water_maxFreq);
end

WV_tc=conj(V_tc - FV_tc);

[WaterSV_tc,~,~] = svd(WV_tc,0);

IOP=diag(ones(NbT,1));

MRSIData_ctkkk =	temp_ckkk_t*(IOP-WaterSV_tc(:,1:FiltParam.Comp)*WaterSV_tc(:,1:FiltParam.Comp)');

MRSIData_ctkkk=permute(reshape(MRSIData_ctkkk,[NbCoil,Nr,Nc,Ns,NbT]),[1,5,2,3,4]);

% figure;
% TEST =  Fast_HSVD_Filter(permute(temp_ckkk_t(1,:),[2,1]),1.0/mrProt.DwellTime,128,FiltParam.Water_minFreq,FiltParam.Water_maxFreq);
% hold on 
% plot(abs(fft(temp_ckkk_t(1,:),[],2)))
% plot(abs(fft(MRSIData_ctkkk(1,:,1,1,1),[],2)))
% plot(abs(fft(TEST)))


end
