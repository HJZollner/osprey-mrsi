function [ReconData_trrr,U_rrrc,V_tc,S,DynamicFreqMap,costFunVal] = LowRankTGV_ModelRecon(mrsiData_ckkkt,mrsiReconParams)


if ~exist([mrsiReconParams.Log_Dir '/LowRankTGV_Recon']);mkdir([mrsiReconParams.Log_Dir '/LowRankTGV_Recon']);end

fprintf('Compute initial SVD Adjoint solution...\n');
SiOri=size(mrsiData_ckkkt);
nDimsOri= ndims(mrsiData_ckkkt);


%apply brain mask for intitial solution
NbT=SiOri(end);
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;
DelayT=mrsiReconParams.AcqDelay;
Time_rrrt=permute(repmat(([0 :(NbT-1)]'/Fs+DelayT),[1 SiOri(2) SiOri(3) SiOri(4)]),[2,3,4,1]);
Freqshift_crrrt= permute(repmat(exp(-2*pi*1i*Time_rrrt.*repmat(mrsiReconParams.WaterFreqMap,[1 1 1 NbT])),[1 1 1 1 SiOri(1)]),[5 1 2 3 4]);
clear  Time_rrrt

temp_ckkkt=ifft(ifft(ifft(mrsiData_ckkkt,[],2),[],3),[],4);% c r r r t
temp_ckkkt=Freqshift_crrrt.*temp_ckkkt;

clear  Freqshift_crrrt
BrainMask_1rrr1=reshape(mrsiReconParams.BrainMask,[1 size(mrsiReconParams.BrainMask) 1 ]);
temp_ckkkt=fft(fft(fft(BrainMask_1rrr1.*temp_ckkkt,[],2),[],3),[],4);
clear  BrainMask_crrrt

[Uorig,Sorig,Vorig] = svd(reshape(temp_ckkkt,[],SiOri(end)),0);
clear temp_ckkkt
V_tc=Vorig(:,1:mrsiReconParams.modelOrder);
S=Sorig(1:mrsiReconParams.modelOrder,1:mrsiReconParams.modelOrder);
VisualizeSpectral( V_tc,S, [mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_Initial']);

U_ckkkc=reshape(Uorig(:,1:mrsiReconParams.modelOrder), SiOri(1),SiOri(2),SiOri(3),SiOri(4),[]);

%% TGV Parameters
kmask=mrsiReconParams.kmask;
alpha = 1E-12*mrsiReconParams.mu_tv; %  
maxit = 1000 ;          % use 1000 Iterations for optimal image quality   
Threshold=1E-7;
Brainmask_1rrr=reshape(mrsiReconParams.BrainMask,[1 size(mrsiReconParams.BrainMask,1) size(mrsiReconParams.BrainMask,2) size(mrsiReconParams.BrainMask,3)]);
reduction = 2^-8;
Init_U_rrrc=zeros(SiOri(2),SiOri(3),SiOri(4),mrsiReconParams.modelOrder);
fprintf('Start initial reconstruction...\n');

parfor k=1:mrsiReconParams.modelOrder
    [Init_U_rrrc(:,:,:,k) e] = tgv2_l2_3D_multiCoil(U_ckkkc(:,:,:,:,k),Brainmask_1rrr.*mrsiReconParams.SENSE, kmask, 2*alpha, alpha, maxit, reduction,Threshold);
end

VisualizeSpectral( V_tc,S, [mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_Initial']);

VisualizeTGV( Init_U_rrrc,[mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_muTV', num2str(mrsiReconParams.mu_tv),'_Initial']);

SENSE_crrrt=repmat(mrsiReconParams.SENSE,[1 1 1 1 SiOri(end)]);
OriginalData_trrr= permute(conj(SENSE_crrrt).*ifft(ifft(ifft(mrsiData_ckkkt,[],2),[],3),[],4),[1 5 2 3 4]);
OriginalData_frrr =fft(squeeze(sum(OriginalData_trrr,1)),[],1);

clear  SENSE_crrrt US_ckkkc U_ckkkc Uorig Vorig OriginalData_trrr IFFTData

U_rrrc=Init_U_rrrc;
SizeVol = size(U_rrrc);

fprintf(makeSectionDisplayHeader('Start the recurrence reconstruction...\n'));

%% TGV Parameters
%kmask=mrsiReconParams.kmask;


maxit = mrsiReconParams.LRTGVModelParams.maxit;        
minit = round(mrsiReconParams.LRTGVModelParams.maxit/3);
Threshold=1;% Not used Anymore

[U_rrrc,V_tc,S, costFunVal,DynamicFreqMap] = tgv2_l2_3D_SENSE_MRSILowRank(mrsiData_ckkkt,U_rrrc,V_tc,S,OriginalData_frrr, 2*mrsiReconParams.mu_tv, mrsiReconParams.mu_tv, maxit,minit,mrsiReconParams,Threshold);

VisualizeTGV(U_rrrc,[mrsiReconParams.Log_Dir,filesep,mrsiReconParams.NameData,'_Final']);
VisualizeSpectral( V_tc,S, [ mrsiReconParams.Log_Dir,filesep,mrsiReconParams.NameData,'_Final'])
s=[ mrsiReconParams.Log_Dir,filesep,mrsiReconParams.NameData,'_Final_FrequencyMap.ps'];
figs=figure('visible','off');imagesc(Vol2Image(DynamicFreqMap),[-50 50]);
colorbar; title('Final Adaptative Frequency Map');print(figs, '-dpsc2',s);

fprintf('Recombination of the Data...\n');
TempRes_rrrt=formTensorProduct(U_rrrc,V_tc*S);
ReconData_trrr=permute(TempRes_rrrt,[ndims( TempRes_rrrt),1:(ndims( TempRes_rrrt)-1)]);

end
