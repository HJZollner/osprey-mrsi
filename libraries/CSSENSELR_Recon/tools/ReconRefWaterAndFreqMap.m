function [ WaterFreqMap_rrr,Water_trrr,Water_rrr]=ReconRefWaterAndFreqMap(mrsiReconParams)



MSize_data=size(mrsiReconParams.mrsiData_ctkkk);
WSize_data=size(mrsiReconParams.Water_ctkkk);

RowSet=[1:ceil(WSize_data(3)/2), ceil(MSize_data(3)-WSize_data(3)/2+1):MSize_data(3)];
ColSet= [1:ceil(WSize_data(4)/2), ceil(MSize_data(4)-WSize_data(4)/2+1):MSize_data(4)];
SlcSet= [1:ceil(WSize_data(5)/2.0), ceil(MSize_data(5)-WSize_data(5)/2.0+1):MSize_data(5)];


% Simple FFT and Coil Combination
reduction = 2^(-8);     % usually there is no need to change this
alpha = 1E-32; %no Regularization
maxit = 1000;            % use 1000 Iterations for optimal image quality

%for k=1:size(US_ckkc,4)
fprintf('Start the Water reconstruction...\n');


SENSE=mrsiReconParams.SENSE;

clear Water_ckkk

	fprintf('Reconstructing Water Signal to MRSI resolution.\n');
	
	[X,Y,Z] = ndgrid(1:MSize_data(3), 1:MSize_data(4),1:MSize_data(5));
	xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;zc=floor(MSize_data(5)/2)+1;
	temp = ((2*(X-xc)/WSize_data(3)).^2 + (2*(Y-yc)/WSize_data(4)).^2+ (2*(Z-zc)/WSize_data(5)).^2).^0.5 ;
	HKernel = fftshift(0.5*(1+cos(pi*temp)));  
	HKernel_11kkk=reshape(HKernel,[1 1 size(HKernel)]);
	
	WaterZeroPad_ctkkk=zeros(WSize_data(1),WSize_data(2), MSize_data(3),MSize_data(4),MSize_data(5));
	WaterZPadKmask=zeros(MSize_data(3),MSize_data(4),MSize_data(5));
	WaterZeroPad_ctkkk(:,:,RowSet,ColSet,SlcSet)=mrsiReconParams.Water_ctkkk;
	WaterZPadKmask(RowSet,ColSet,SlcSet)=mrsiReconParams.Water_kmask;
	
	WaterZeroPad_ctkkk  = single(WaterZeroPad_ctkkk); %.*HKernel_11kkk;

	[Uorig,Sorig,Vorig] = svd(reshape(permute(WaterZeroPad_ctkkk,[1 3 4 5 2]),[],WSize_data(2)),0);

	SingVal=diag(Sorig)./cumsum(diag(Sorig));%(abs(Sorig(1,1));
	OrderWater=min(find(abs(diff(SingVal))<1E-3));

	clear HKernel_ckkk
	if OrderWater>mrsiReconParams.modelOrder
	    OrderWater=mrsiReconParams.modelOrder;
	end
	if OrderWater<8
	    OrderWater=8;
	end
	V=Vorig(:,1:OrderWater);
	U_ckkkc=reshape(Uorig, WSize_data(1),MSize_data(3),MSize_data(4),MSize_data(5),[]);
	Recon_US_rrrc=zeros(MSize_data(3),MSize_data(4),MSize_data(5),OrderWater);

	%alpha = 5E-6;Threshold=1E-6;
	
	
	parfor k=1:OrderWater
		alpha = 5E-6*norm( vectorizeArray(U_ckkkc(:,:,:,:,k)));Threshold=1E-9;maxit=500;
		[Recon_US_rrrc(:,:,:,k) e] = tgv2_l2_3D_multiCoil(U_ckkkc(:,:,:,:,k),SENSE, WaterZPadKmask, 2*alpha, alpha, maxit, reduction,Threshold);
		Recon_US_rrrc(:,:,:,k)=Sorig(k,k)*Recon_US_rrrc(:,:,:,k);
	end

	Water_trrr=formTensorProduct(Recon_US_rrrc,V);
	Water_trrr=permute(Water_trrr,[4,1,2,3]);



%Water_rrr=squeeze(mean(Water_trrr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),1));
alpha = 5E-6*norm(WaterZeroPad_ctkkk(:));Threshold=1E-9;maxit=500;
[Water_rrr e] = tgv2_l2_3D_multiCoil(squeeze(sum(WaterZeroPad_ctkkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),2)),SENSE, WaterZPadKmask, 2*alpha, alpha, maxit, reduction,Threshold);


if ~mrsiReconParams.SetInitFreqMapToZero
    if mrsiReconParams.LoadNiftiFreqMap 
    	[WaterFreqMap_rrr ~]=BuildFieldmap('./B0MapMagn.nii.gz','./B0MapPhase.nii.gz', 2.080, 3.100);
        if sum(size(WaterFreqMap_rrr)==MSize_data(3:5))<3
             size_WaterFreqMap_rrr = size(WaterFreqMap_rrr)
             size_MRSI_data = MSize_data(3:5)
             error('B0 field map loaded from the nifti files doesnt have the same dimensions as the MRSI data!')
	end
    else
	WaterFreqMap_rrr=DetermineSingleLowFreq( Water_trrr,60,mrsiReconParams);
    end
else
    WaterFreqMap_rrr=0*Water_rrr;
end

WaterFreqMap_rrr=WaterFreqMap_rrr-mean(WaterFreqMap_rrr(mrsiReconParams.BrainMask>0));%remove intercept
close all;figs=figure('visible', 'off');

s=[mrsiReconParams.Log_Dir,filesep,mrsiReconParams.NameData, '_WB_Water_Amp_Phase_FreqMap.ps'];
if exist(s);delete(s);end
figs=figure('visible', 'off');

 plotImage= Vol2Image( flip(permute( Water_rrr,[2,1,3]),1) );

 imagesc(abs(plotImage),[0 max(abs(Water_rrr(:)))]);daspect([1 1 1]);
colormap default;colorbar;title('Water Amp Recon (Used for Normalization)')
print(figs, '-append', '-dpsc2', s);

Water2_rrr=squeeze(mean(Water_trrr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),1));
 plotImage= Vol2Image( flip(permute( Water2_rrr ,[2,1,3]) ,1) );
 imagesc(abs(plotImage),[0 max(abs(Water2_rrr(:)))]);daspect([1 1 1]);

colormap default;colorbar;title('Water Signal Recon Amplitude (Should be consitent with Water Amp Recon)')
print(figs, '-append', '-dpsc2', s);

imagesc(angle(plotImage));daspect([1 1 1]);
colormap default;colorbar;title('Water Signal Recon Phase')
print(figs, '-append', '-dpsc2', s);

 plotImage= Vol2Image( flip(permute(WaterFreqMap_rrr,[2,1,3]),1) );
imagesc(plotImage,[-50 50]);daspect([1 1 1]);
colormap default;colorbar;title('Freq Map Based on Water Signal')
print(figs, '-append', '-dpsc2', s);

 plotImage= Vol2Image( flip(permute(mrsiReconParams.BrainMask.*WaterFreqMap_rrr,[2,1,3]) ,1) ) ;
imagesc(plotImage,[-50 50]);daspect([1 1 1]);
colormap default;colorbar;title('Brain-masked Freq Map Based on Water Signal')
print(figs, '-append', '-dpsc2', s);


close all;


fprintf('Coil profile computation and water reconstruction finished.\n');
end
