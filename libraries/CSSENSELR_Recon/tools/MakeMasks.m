function [ ImMask, BrainMask,SkMask ] = MakeMasks(mrsiReconParams )

FiltParam.FilterFreq=0;%5;%hz
FiltParam.Water_minFreq=-60;%-60;%hz
FiltParam.Water_maxFreq=60;%60hz
FiltParam.Comp=16;

mrProt=mrsiReconParams.Water_mrProt;
kmask=mrsiReconParams.Water_kmask;

MSize_data=size(mrsiReconParams.mrsiData_ctkkk);
WSize_data=size(mrsiReconParams.Water_ctkkk);

RowSet=[1:ceil(WSize_data(3)/2), ceil(MSize_data(3)-WSize_data(3)/2+1):MSize_data(3)];
ColSet= [1:ceil(WSize_data(4)/2), ceil(MSize_data(4)-WSize_data(4)/2+1):MSize_data(4)];
SlcSet= [1:ceil(WSize_data(5)/2.0), ceil(MSize_data(5)-WSize_data(5)/2.0+1):MSize_data(5)];

Water_ckkk=sqz(mean(mrsiReconParams.Water_ctkkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),2));


[X,Y,Z] = ndgrid(1:MSize_data(3), 1:MSize_data(4),1:MSize_data(5));
xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;zc=floor(MSize_data(5)/2)+1;
temp = ((2*(X-xc)/WSize_data(3)).^2 + (2*(Y-yc)/WSize_data(4)).^2+ (2*(Z-zc)/WSize_data(5)).^2).^0.5 ;
HKernel = fftshift(0.5*(1+cos(pi*temp)));  

	
WaterZeroPad_ckkk=zeros( WSize_data(1),MSize_data(3),MSize_data(4),MSize_data(5));
WaterZPadKmask=zeros(MSize_data(3),MSize_data(4),MSize_data(5));
WaterZeroPad_ckkk(:,RowSet,ColSet,SlcSet)=Water_ckkk;
WaterZPadKmask(RowSet,ColSet,SlcSet)=mrsiReconParams.Water_kmask;
WaterZeroPad_ckkk  = single(WaterZeroPad_ckkk);%.*HKernel;

WaterZeroPad_rrr= sqz(sum(abs(ifft(ifft(ifft(WaterZeroPad_ckkk,[],2),[],3),[],4)).^2,1));

clear WaterZeroPad_ckkk

MRSI_frrr=sqz(sum( abs(fft(ifft(ifft(ifft(mrsiReconParams.mrsiData_ctkkk,[],3),[],4),[],5),[],2)).^2 ,1) );



[ ~ , low_bnd_L ]=min(abs( mrsiReconParams.LipidMinPPM - mrsiReconParams.ppm));
[ ~ , high_bnd_L ]=min(abs( mrsiReconParams.LipidMaxPPM - mrsiReconParams.ppm));

Lipid_rrr=squeeze(sum(abs(MRSI_frrr(low_bnd_L:high_bnd_L,:,:,:)),1));
Lipid_rrr=Lipid_rrr/max(Lipid_rrr(:));
ImSize=MSize_data(3:5);
CornerMask=zeros(ImSize);
RowMask=[1:round(ImSize(1)/10),(ImSize(1)-round(ImSize(1)/10)+1):ImSize(1)];
ColMask=[1:round(ImSize(2)/10),(ImSize(2)-round(ImSize(2)/10)+1):ImSize(2)];
CornerMask(RowMask,ColMask,1:end)=1;
CenterMask=zeros(ImSize);
RowMask=[round(ImSize(1)*0.4):round(ImSize(1)*0.6)];
ColMask=[round(ImSize(2)*0.4):round(ImSize(2)*0.6)];
CenterMask(RowMask,ColMask,2:end-1)=1;

WaterThres=0.5*(quantile(WaterZeroPad_rrr(CornerMask(:)>0),0.95)+quantile(WaterZeroPad_rrr(CenterMask(:)>0),0.25));
ImMask=WaterZeroPad_rrr>WaterThres;
%1st estimate
LipidThres=0.5*(quantile(Lipid_rrr(CornerMask(:)>0),0.95)+0.5*quantile(Lipid_rrr(:),mrsiReconParams.AutoMask_Quantile));
%estimation refinement
LipidThres=0.5*(quantile(Lipid_rrr(CornerMask(:)>0),0.95)+quantile(Lipid_rrr(Lipid_rrr(:)>LipidThres),0.25));
%SkMask=Lipid_rrr>LipidThres;


for RedExp=1:20;
    SkMask=Lipid_rrr>LipidThres*(1/2^(RedExp-10));

    CC= bwconncomp(SkMask); 
   NbComp(RedExp) = CC.NumObjects;
   NbVox(RedExp) = sum(SkMask(:)>0);
end
[~, MaxNbComp] = max(NbComp);
[~, MaxCurvI] = min(NbComp(1:MaxNbComp)./NbVox(1:MaxNbComp));

RedExp=MaxCurvI;
SkMask=Lipid_rrr>LipidThres*(1/2^(RedExp-10));

s=[mrsiReconParams.Log_Dir,filesep,mrsiReconParams.NameData,'_AutoMasks_Lipid_Head.ps'];
if exist(s);delete(s);end
figs=figure('visible','off');
volimagesc(WaterZeroPad_rrr);%,[ 0, 10*mean(image2plot(:))] );
title('Water intensity map');
print(figs, '-append', '-dpsc2', s);

volimagesc(Lipid_rrr);%,[ 0, 10*mean(image2plot(:))] );
title('Lipids intensity map');
print(figs, '-append', '-dpsc2', s);

volimagesc(ImMask);
title('1st image mask estimate');
print(figs, '-append', '-dpsc2', s);
volimagesc(SkMask);
title('1st lipid mask estimate');
print(figs, '-append', '-dpsc2', s); 

seD = strel('sphere',1);
seD_2D = strel('disk',1);
seC = strel('sphere',round(max(size(SkMask))/2));
seC_2D = strel('disk',round(max(size(SkMask))/2));
seO5 = strel('sphere',5);
seO5_2D = strel('disk',5);

%Padding the Masks
Temp(:,:,1)=ImMask(:,:,1);
Temp(:,:,size(ImMask,3)+2)=ImMask(:,:,end);
Temp(:,:,2:end-1)=ImMask;
ImMask=Temp;

%Padding the Masks
Temp(:,:,1)=SkMask(:,:,1);
Temp(:,:,size(SkMask,3)+2)=SkMask(:,:,end);
Temp(:,:,2:end-1)=SkMask;
SkMask=Temp;

BrainMask=0*ImMask;
DilSkMask=0*SkMask;

ImMask = imopen(imclose(ImMask,seO5),seO5);
%closing the skull mask in 3D is not a good idea because it is open n the bottom
for z=1:size(SkMask,3)
    DilSkMask(:,:,z) = imclose(SkMask(:,:,z),seC_2D);
end
BrainMask =  imopen(imclose(DilSkMask -SkMask,seO5),seO5);

ImMask = (ImMask| BrainMask|SkMask );
ImMask= imopen(imclose(ImMask,seO5),seO5);

SkMask=single(SkMask(:,:,2:end-1));
BrainMask=single(BrainMask(:,:,2:end-1));
ImMask=single(ImMask(:,:,2:end-1));

volimagesc(SkMask);
title('Lipid mask estimate');
print(figs, '-append', '-dpsc2', s); 
volimagesc(ImMask);
title('Image mask estimate');
print(figs, '-append', '-dpsc2', s); 
volimagesc(BrainMask);
title('Brain mask estimate');
print(figs, '-append', '-dpsc2', s); 
end

