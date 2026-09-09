function [SENSE, WaterFreqMap_rrr,Water_trrr,Water_rrr,UnCorrSENSE]=ComputeSENSEProfilesWithESPIRiT(mrsiReconParams)


ksize = mrsiReconParams.ESPIRIT_kernel;

% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.
eigThresh_1 = 0.02;
% threshold of eigen vector decomposition in image space.
eigThresh_2 = 0.95;

MSize_data=size(mrsiReconParams.mrsiData_ctkkk);
WSize_data=size(mrsiReconParams.Water_ctkkk);

RowSet=[1:ceil(WSize_data(3)/2), ceil(MSize_data(3)-WSize_data(3)/2+1):MSize_data(3)];
ColSet= [1:ceil(WSize_data(4)/2), ceil(MSize_data(4)-WSize_data(4)/2+1):MSize_data(4)];
SlcSet= [1:ceil(WSize_data(5)/2.0), ceil(MSize_data(5)-WSize_data(5)/2.0+1):MSize_data(5)];

DATA=permute(sum(mrsiReconParams.Water_ctkkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),2),[3,4,5,1,2]);
DATA=fft(fft(fft(fftshift(fftshift(fftshift(ifft(ifft(ifft(DATA,[],1),[],2),[],3),1),2),3),[],1),[],2),[],3);
DATA=fftshift(fftshift(fftshift(DATA,1),2),3);
Check_kmask=fftshift(fftshift(fftshift(mrsiReconParams.Water_kmask,1),2),3);

[Nc,~,sx,sy,sz] = size(mrsiReconParams.mrsiData_ctkkk);

%% Compute ESPIRiT EigenVectors

% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels

count=0;A=[];

for z=1:(size(DATA,3)-ksize(3)+1)
    for y=1:(size(DATA,2)-ksize(2)+1)
        for x=1:(size(DATA,1)-ksize(1)+1)
            Check_Kernel=Check_kmask(x:(x+ksize(1)-1),y:(y+ksize(2)-1),z:(z+ksize(3)-1) ) ;
            if ( mean(Check_Kernel(:))==1) %then kernel is acquired at this location
                count = count+1;

                A(count,:,:) = reshape(DATA(x:(x+ksize(1)-1),y:(y+ksize(2)-1),z:(z+ksize(3)-1),:),1,ksize(1)*ksize(2)*ksize(3),Nc);
            end
        end
    end
end

fprintf([num2str(count), ' kernels of size ',num2str(ksize(1)),'x',num2str(ksize(2)),'x',num2str(ksize(3)) , ' found in the Water measurement data.\n']);
filling=1;
while count < 50
    fprintf(['Too few number kernels. Computing again with a ',num2str(filling),' filling threshold.\n']);
    count=0;A=[];filling=filling*0.95;

    for z=1:(size(DATA,3)-ksize(3)+1)
        for y=1:(size(DATA,2)-ksize(2)+1)
    		for x=1:(size(DATA,1)-ksize(1)+1)
    		    Check_Kernel=Check_kmask(x:(x+ksize(1)-1),y:(y+ksize(2)-1),z:(z+ksize(3)-1) ) ;
    		    if ( mean(Check_Kernel(:))>=filling) %then kernel is acquired at this location
    		        count = count+1;

    		        A(count,:,:) = reshape(DATA(x:(x+ksize(1)-1),y:(y+ksize(2)-1),z:(z+ksize(3)-1),:),1,ksize(1)*ksize(2)*ksize(3),Nc);
    		    end
    		end
        end
    end

    fprintf([num2str(count), ' kernels of size ',num2str(ksize(1)),'x',num2str(ksize(2)),'x',num2str(ksize(3)) , ' found in the Water measurement data with ',num2str(filling),' kernel filling.\n']);
end
A = reshape(A,size(A,1),size(A,2)*size(A,3));
[~,S,V] = svd(A);

k = reshape(V,ksize(1),ksize(2),ksize(3),Nc,size(V,2));
S = diag(S);S = S(:);
clear V A

idx = max(find(S >= S(1)*eigThresh_1));


%%
% crop kernels and compute eigen-value decomposition in image space to get
% maps
[M,W] = kernelEig(k(:,:,:,:,1:idx),[sx,sy,sz]);

%%
% crop sensitivity maps
%
if (mean(mrsiReconParams.ImMask(:))==0 || mean(mrsiReconParams.ImMask(:))==1 ) % then ImMask is useless and we use the weight to mask the coil sensitivity
    maps = M(:,:,:,:,end).*repmat(W(:,:,:,end)>eigThresh_2,[1,1,1,Nc]);
else
    maps = M(:,:,:,:,end).*repmat(mrsiReconParams.ImMask,[1,1,1,Nc]);
end

SENSE=permute(maps,[4,1,2,3]);

Data_name=mrsiReconParams.NameData;
s1=[mrsiReconParams.Log_Dir,filesep,'ESPIRIT_CoilProfiles_FromWaterMeas_Abs_',Data_name,'.ps'];
s2=[mrsiReconParams.Log_Dir,filesep,'ESPIRIT_CoilProfiles_FromWaterMeas_Phase_',Data_name,'.ps'];
s5=[mrsiReconParams.Log_Dir,filesep,'Water_Homogenization_',Data_name,'.ps'];
if exist(s1);delete(s1);end
if exist(s2);delete(s2);end
if exist(s5);delete(s5);end

figs=figure('visible', 'off');


SENSE(find(isnan(SENSE)))=0;



reduction = 2^(-8);     % usually there is no need to change this
fprintf('Start the Water reconstruction...\n');

UnCorrSENSE=SENSE;

fprintf('Reconstructing Water Signal to MRSI resolution.\n');


WaterZeroPad_ctkkk=zeros(WSize_data(1),WSize_data(2), MSize_data(3),MSize_data(4),MSize_data(5));
WaterZPadKmask=zeros(MSize_data(3),MSize_data(4),MSize_data(5));
WaterZeroPad_ctkkk(:,:,RowSet,ColSet,SlcSet)=mrsiReconParams.Water_ctkkk;
WaterZPadKmask(RowSet,ColSet,SlcSet)=mrsiReconParams.Water_kmask;

WaterZeroPad_ctkkk  = single(WaterZeroPad_ctkkk);

[Uorig,Sorig,Vorig] = svd(reshape(permute(WaterZeroPad_ctkkk,[1 3 4 5 2]),[],WSize_data(2)),0);

SingVal=diag(Sorig)./cumsum(diag(Sorig));
OrderWater=min(find(abs(diff(SingVal))<1E-3));

if OrderWater>mrsiReconParams.modelOrder
    OrderWater=mrsiReconParams.modelOrder;
end
if OrderWater<8
    OrderWater=8;
end
V=Vorig(:,1:OrderWater);
U_ckkkc=reshape(Uorig, WSize_data(1),MSize_data(3),MSize_data(4),MSize_data(5),[]);
Recon_US_rrrc=zeros(MSize_data(3),MSize_data(4),MSize_data(5),OrderWater);

parfor k=1:OrderWater
    alpha = 5E-6*norm( vectorizeArray(U_ckkkc(:,:,:,:,k)));Threshold=1E-9;maxit=500;
    [Recon_US_rrrc(:,:,:,k) ~] = tgv2_l2_3D_multiCoil(U_ckkkc(:,:,:,:,k),SENSE, WaterZPadKmask, 2*alpha, alpha, maxit, reduction,Threshold);
    Recon_US_rrrc(:,:,:,k) = Sorig(k,k)*Recon_US_rrrc(:,:,:,k);
end

Water_trrr=formTensorProduct(Recon_US_rrrc,V);
Water_trrr=permute(Water_trrr,[4,1,2,3]);

alpha = 5E-6*norm(WaterZeroPad_ctkkk(:));Threshold=1E-9;maxit=500;
[Water_rrr ~] = tgv2_l2_3D_multiCoil(squeeze(sum(WaterZeroPad_ctkkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),2)),SENSE, WaterZPadKmask, 2*alpha, alpha, maxit, reduction,Threshold);


if ~mrsiReconParams.SetInitFreqMapToZero
    if mrsiReconParams.LoadNiftiFreqMap
        [WaterFreqMap_rrr ~]=BuildFieldmap('./B0MapMagn.nii.gz','./B0MapPhase.nii.gz', 2.080, 3.100);
        if sum(size(WaterFreqMap_rrr)==MSize_data(3:5))<3
            size_WaterFreqMap_rrr = size(WaterFreqMap_rrr)
            size_MRSI_data = MSize_data(3:5)
            error('B0 field map loaded from the nifti files doesnt have the same dimensions as the MRSI data!')
    	end
    else
        WaterFreqMap_rrr=DetermineSingleLowFreq( Water_trrr,mrsiReconParams.WaterRefFreqSearchRange,mrsiReconParams);
    end
else
    WaterFreqMap_rrr=0*Water_rrr;
end


close all;figs=figure('visible', 'off');
for k=1:size(mrsiReconParams.Water_ctkkk,1)

    plotImage= Vol2Image( squeeze(SENSE(k,:,:,:)) );

    imagesc(abs(squeeze(plotImage)));
    colormap default;
    print(figs, '-append', '-dpsc2', s1);

    imagesc(angle(plotImage ));
    colormap default;
    print(figs, '-append', '-dpsc2', s2);

end

s=[mrsiReconParams.Log_Dir,filesep,mrsiReconParams.NameData, '_Water_Amp_Phase_FreqMap.ps'];
if exist(s);delete(s);end
figs=figure('visible', 'off');

plotImage= Vol2Image( flip(permute( Water_rrr,[2,1,3]),1) );

imagesc(abs(plotImage),[0 max(abs(Water_rrr(:)))]);daspect([1 1 1]);
colormap default;colorbar;title('Water Amp Recon')
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


close all;


fprintf('Coil profile computation and water reconstruction finished.\n');
end
