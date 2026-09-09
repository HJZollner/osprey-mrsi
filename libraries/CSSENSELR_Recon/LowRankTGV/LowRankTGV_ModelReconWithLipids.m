function [ReconDataShort_trrr,ReconMetab_trrr,ReconLipid_trrr,U_rrrc,V_tc,S,LipidRM_ff,DynamicFreqMap] = LowRankTGV_ModelReconWithLipids(mrsiData_ckkkt,RecParams)

if ~exist([RecParams.Log_Dir '/LowRankTGV_Recon']);mkdir([RecParams.Log_Dir '/LowRankTGV_Recon']);end

fprintf(makeSectionDisplayHeader('Start the iterative reconstruction...\n'));

%% TGV Parameters
alpha = RecParams.mu_tv;
maxit = RecParams.LRTGVModelParams.maxit;        
minit = round(RecParams.LRTGVModelParams.maxit/3);

Threshold=1E-9;

[U_rrrc,V_tc,S,ReconLipid_rrrt,MetabData_rrrt,DynamicFreqMap,costFunVal,LipidRM_ff] = tgv2_l2_3D_SENSE_LowRank_MetabLipCoRecon(mrsiData_ckkkt, 2*alpha, alpha, maxit,minit,RecParams,Threshold);

VisualizeTGV( U_rrrc,[RecParams.Log_Dir,filesep,RecParams.NameData,'_Final']);
VisualizeSpectral( V_tc,S, [ RecParams.Log_Dir,filesep,RecParams.NameData,'_Final'])
s=[ RecParams.Log_Dir,filesep,RecParams.NameData,'_Final_FrequencyMap.ps'];
figs=figure('visible','off');volimagesc(DynamicFreqMap);
colorbar; title('Final Adaptative Frequency Map');print(figs, '-dpsc2',s);

fprintf('Data Recombination ...\n');
TempRes_rrrt=formTensorProduct(U_rrrc,V_tc*S);
ReconDataShort_trrr=permute(TempRes_rrrt,[ndims( TempRes_rrrt),1:(ndims( TempRes_rrrt)-1)]);

ReconMetab_frrr = permute(fft(MetabData_rrrt,[],4),[4,1,2,3]);
ReconMetab_trrr=ifft(ReconMetab_frrr,[],1);

ReconLipid_trrr = permute(ReconLipid_rrrt,[ndims( ReconLipid_rrrt),1:(ndims( ReconLipid_rrrt)-1)]);

end
