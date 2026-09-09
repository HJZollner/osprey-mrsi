function [SENSE, WaterFreqMap_rrr,Water_trrr,Water_rrr,UnCorrSENSE]=ComputeSENSEProfilesFromWaterAcq(mrsiReconParams)

%Water_crr=squeeze(sum(ifft(ifft(mrsiReconParams.Water_ctkk(:,1:round(end/10),:,:),[],3),[],4),2));

MSize_data=size(mrsiReconParams.mrsiData);
WSize_data=size(mrsiReconParams.Water_ctkkk);

nDimsOri= ndims(mrsiReconParams.Water_ctkkk);

%}
% *************************************************************************
% Compute Profile in Original Resolution
% *************************************************************************

Water_ctrrr = ifft(ifft(ifft(mrsiReconParams.Water_ctkkk,[],3),[],4),[],5);

[Uorig,Sorig,Vorig] = svd(reshape(permute(Water_ctrrr,[1 3 4 5 2]),[],WSize_data(2)),0);
USHead_crrrc=reshape(Uorig*Sorig, WSize_data(1),WSize_data(3),WSize_data(4),WSize_data(5),[]);
WaterAmpOrig_crrr=abs(USHead_crrrc(:,:,:,:,1));
WaterPhOrig_crrr=angle(USHead_crrrc(:,:,:,:,1));



% *************************************************************************
% Compute Profile in High Resolution
% *************************************************************************

RowSet=[1:WSize_data(3)/2, (MSize_data(3)-WSize_data(3)/2+1):MSize_data(3)];
ColSet= [1:WSize_data(4)/2, (MSize_data(4)-WSize_data(4)/2+1):MSize_data(4)];
SlcSet= [1:WSize_data(5)/2, (MSize_data(5)-WSize_data(5)/2+1):MSize_data(5)];

%[Uorig,Sorig,Vorig] = svd(reshape(permute(BWater_ctrr,[1 3 4 2]),[],BWSize_data(2)),0);
%USBody_crrc=reshape(Uorig*Sorig, BWSize_data(1),BWSize_data(3),BWSize_data(4),[]);

%Gaussian Kernel
% sigma=(size(mrsiReconParams.Water_ctkkk,3)+size(mrsiReconParams.Water_ctkkk,4))/6;%3; %Optimized with comparison to the Full Res Water Data
% [X,Y,Z] = ndgrid(1:MSize_data(3), 1:MSize_data(4),1:MSize_data(5));
% xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;zc=floor(MSize_data(5)/2)+1;
% exponent = -((X-xc).^2 + (Y-yc).^2+ (Z-zc).^2)./(2*sigma^2);
%Kernel = fftshift(exp(exponent)); % no need to normalize

%Hanning Kernel
sigma=(size(mrsiReconParams.Water_ctkkk,3)+size(mrsiReconParams.Water_ctkkk,4))/6;%3; %Optimized with comparison to the Full Res Water Data
[X,Y,Z] = ndgrid(1:MSize_data(3), 1:MSize_data(4),1:MSize_data(5));
xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;zc=floor(MSize_data(5)/2)+1;
temp = ((2*(X-xc)/WSize_data(3)).^2 + (2*(Y-yc)/WSize_data(4)).^2+ (2*(Z-zc)/WSize_data(5)).^2).^0.5 ;
HKernel = fftshift(0.5*(1+cos(pi*temp)));  

% Zero Padding Water and filtering
WaterZeroPad_ctkkk=zeros(WSize_data(1),WSize_data(2), MSize_data(3),MSize_data(4),MSize_data(5));
WaterZeroPad_ctkkk(:,:,RowSet,ColSet,SlcSet)=mrsiReconParams.Water_ctkkk;
HKernel_ctkkk=permute(repmat(HKernel,[1,1,1,WSize_data(1),WSize_data(2)]),[4 5 1 2 3]);
WaterZeroPad_ctrrr  = ifft(ifft(ifft(WaterZeroPad_ctkkk.*HKernel_ctkkk,[],3),[],4),[],5);

WaterAmp_crrr=squeeze(mean(abs(WaterZeroPad_ctrrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:)),2));
WaterPh_crrr=squeeze(angle(sum(WaterZeroPad_ctrrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),2)));
Water_crrr=WaterAmp_crrr.*exp(1j*WaterPh_crrr);
%BodyNorm=norm(Water_crrr(:));
SENSESQ = squeeze(sqrt(sum(abs(Water_crrr(:,:,:,:)).^2,1)));
SENSE = Water_crrr;
clear Water_ctrrr  WaterZeroPad_ctrrr HKernel_ctkkk
%%

Data_name=mrsiReconParams.NameData;
s1=[mrsiReconParams.Log_Dir,filesep,'WHB_SENSECoilProfiles_FromWaterMeas_Abs_',Data_name,'.ps'];
s2=[mrsiReconParams.Log_Dir,filesep,'WHB_SENSECoilProfiles_FromWaterMeas_Phase_',Data_name,'.ps'];
s3=[mrsiReconParams.Log_Dir,filesep,'WHB_HeadWater_perCoil_Abs_',Data_name,'.ps'];
s4=[mrsiReconParams.Log_Dir,filesep,'WHB_HeadWater_perCoil_Phase_',Data_name,'.ps'];
s5=[mrsiReconParams.Log_Dir,filesep,'Water_Homogenization_',Data_name,'.ps'];
if exist(s1);delete(s1);end
if exist(s2);delete(s2);end
if exist(s3);delete(s3);end
if exist(s4);delete(s4);end
if exist(s5);delete(s5);end

figs=figure('visible', 'off');


%SignalThres=mean(SENSESQ(mrsiReconParams.ImMask2D==0))*2;
%SignalThres=quantile(SENSESQ(mrsiReconParams.ImMask==0),0.95)*2;
%SignalThres=mean(SENSESQ(:))*0.05;
%if isnan(SignalThres)
    SignalThres=quantile(SENSESQ(:),1/3);
%end
for k=1:WSize_data(1) 
    
   %SENSE(k,:,:)=ifft(ifft(Kernel.*squeeze(fft(fft(SENSE(k,:,:),[],2),[],3)),[],1),[],2);
    
    SENSE(k,:,:,:)=squeeze(SENSE(k,:,:,:))./SENSESQ;%.*(SENSESQ>SignalThres);
    %SENSE(k,:,:,:)=squeeze(SENSE(k,:,:,:)).*mrsiReconParams.BrainMask;%+~mrsiReconParams.BrainMask*1E-6);
    SENSE(k,:,:,:)=squeeze(SENSE(k,:,:,:)).*mrsiReconParams.ImMask;
    
end
SENSE(find(isnan(SENSE)))=0;
clear Water_crrr exponent
UnCorrSENSE=SENSE;

% Simple FFT and Coil Combination
reduction = 2^(-8);     % usually there is no need to change this
alpha = 1E-32; %no Regularization
maxit = 1000;            % use 1000 Iterations for optimal image quality

%for k=1:size(US_ckkc,4)
fprintf('Start the Water reconstruction...\n');

%Gaussian Kernel
%NX=MSize_data(3);NY= MSize_data(4);NZ= MSize_data(5);
%sigma=mrsiReconParams.GaussianSigma;%1 was not enough
%[X,Y,Z] = ndgrid(1:NX, 1:NY,1:NZ);
%xc=floor(NX/2)+1;yc=floor(NY/2)+1;zc=floor(NZ/2)+1;
%exponent = -((X-xc).^2 + (Y-yc).^2+ (Z-zc).^2)./(2*sigma^2);
%Kernel = fftshift(fftshift(fftshift(exp(exponent)/sum(exp(exponent(:))),1),2),3); % no need to normalize
%Kernel=fft(fft(fft(Kernel,[],1),[],2),[],3);

%smoothing of Maps
NX=MSize_data(3);NY=MSize_data(4);
sigma=mrsiReconParams.GaussianSigma/2;%3 was not enough
[X,Y] = ndgrid(1:NX, 1:NY);
xc=floor(NX/2)+1;yc=floor(NY/2)+1;
exponent = -((X-xc).^2 + (Y-yc).^2)./(2*sigma^2);

Kernel = fftshift(fftshift(exp(exponent)/sum(exp(exponent(:))),1),2); % no need to normalize
Kernel=fft(fft(Kernel,[],1),[],2);

BrainMask_crrr=permute(repmat(mrsiReconParams.BrainMask,[1 1 1 size(SENSE,1)]),[4,1,2,3]);

HKernel_ckkk=permute(repmat(HKernel,[1,1,1,WSize_data(1)]),[4 1 2 3]);
WaterAmp_ckkk=squeeze(mean(abs(WaterZeroPad_ctkkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:)),2));
WaterPh_ckkk=squeeze(angle(sum(WaterZeroPad_ctkkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),2)));
Water_ckkk=WaterAmp_ckkk.*exp(1j*WaterPh_ckkk).*HKernel_ckkk;
FullkSpace=ones(MSize_data(3),MSize_data(4),MSize_data(5));
RealKmask = log(squeeze(sum(abs(Water_ckkk+eps),1)));
RealKmask=RealKmask>(max(RealKmask(:))-15);
clear WaterAmp_ckkk WaterPh_ckkk 
Threshold=1E-8;
WaterHomo=1;
WH_it=1;


RadiusErode=1;
for slc=1:size(mrsiReconParams.BrainMask,3)
    MaskWaterHom(:,:,slc)=1+imerode(mrsiReconParams.BrainMask(:,:,slc),offsetstrel(ones(2*RadiusErode+1)));
end

MaskWaterAdjust=mrsiReconParams.BrainMask;

SmoWater_rrr=zeros(MSize_data(3),MSize_data(4),MSize_data(5));
Water_rrr=zeros(MSize_data(3),MSize_data(4),MSize_data(5));
while WaterHomo>1E-3 & WH_it<10;
    [Water_rrr e] = tgv2_l2_3D_multiCoil(Water_ckkk,SENSE, FullkSpace, 2*alpha, alpha, maxit, reduction,Threshold);

    MeanWater=mean(abs(Water_rrr(MaskWaterAdjust>0)));
    Water_rrr=MaskWaterAdjust.*abs(Water_rrr)/MeanWater+(~MaskWaterAdjust);%.*abs(Water_rrr);
    
    SENSE=SENSE.*BrainMask_crrr;

    for s=1:size(Water_rrr,3)
       % SmoWater_rrr(:,:,s)=abs(ifft(ifft(Kernel.*fft(fft(squeeze(Water_rrr(:,:,s)),[],1),[],2),[],1),[],2));
    end
     
    SmoWater_rrr=Water_rrr;
    SmoWater_rrr=mrsiReconParams.BrainMask.*SmoWater_rrr+(~mrsiReconParams.BrainMask);
    
    subplot(1,2,1);
    imagesc( Vol2Image( Water_rrr ));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;title(['Water Signal It:',num2str(WH_it)]);
    subplot(1,2,2);
    imagesc( Vol2Image(  SmoWater_rrr ));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;title(['Water Signal Smoothed It:',num2str(WH_it)]);
    print(figs,'-bestfit', '-append', '-dpsc2', s5);

    for k=1:WSize_data(1)
        SENSE(k,:,:,:)=squeeze(SENSE(k,:,:,:)).*SmoWater_rrr; %No correction for Synthetic Data 
    end
    WaterHomo=std(SmoWater_rrr(MaskWaterHom>0))/mean(SmoWater_rrr(MaskWaterHom>0)) %0 for Synthetic Data 
    WH_it=WH_it+1
end

clear Water_ckkk

[Uorig,Sorig,Vorig] = svd(reshape(permute(WaterZeroPad_ctkkk,[1 3 4 5 2]).*HKernel_ckkk,[],WSize_data(2)),0);
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


parfor k=1:OrderWater
        [Recon_US_rrrc(:,:,:,k) e] = tgv2_l2_3D_multiCoil(U_ckkkc(:,:,:,:,k),SENSE, FullkSpace, 2*alpha, alpha, maxit, reduction,Threshold);
        Recon_US_rrrc(:,:,:,k)=Sorig(k,k)*Recon_US_rrrc(:,:,:,k);%.*(mrsiReconParams.BrainMask + (~mrsiReconParams.BrainMask)*1E6);
end

Water_trrr=formTensorProduct(Recon_US_rrrc,V,2);
Water_trrr=permute(Water_trrr,[4,1,2,3]);
    
WaterAmp_crrr=squeeze(mean(abs(Water_trrr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:)),1));
WaterPh_crrr=squeeze(angle(sum(Water_trrr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),1)));
Water_rrr=WaterAmp_crrr.*exp(1j*WaterPh_crrr);

WaterFreqMap_rrr=DetermineSingleLowFreq( Water_trrr,60,mrsiReconParams);
WaterFreqMap_rrr=WaterFreqMap_rrr-mean(WaterFreqMap_rrr(mrsiReconParams.BrainMask>0));%remove intercept
close all;figs=figure('visible', 'off');
for k=1:size(mrsiReconParams.Water_ctkkk,1)
    
    plotImage= Vol2Image( squeeze(SENSE(k,:,:,:)) );
    
    imagesc(abs(squeeze(plotImage)));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;
    print(figs, '-append', '-dpsc2', s1);
    
    imagesc(angle(plotImage ));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;
    print(figs, '-append', '-dpsc2', s2);
    
     plotImage= Vol2Image( squeeze(WaterAmpOrig_crrr(k,:,:,:)) );
      imagesc(abs(squeeze(plotImage)));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;
    print(figs, '-append', '-dpsc2', s3);
    
     plotImage= Vol2Image( squeeze(WaterPhOrig_crrr(k,:,:,:)) );
    imagesc(plotImage);%,[ 0, 10*mean(image2plot(:))] )
    colormap default;
    print(figs, '-append', '-dpsc2', s4);    
end

s=[mrsiReconParams.Log_Dir,filesep,mrsiReconParams.NameData, '_WB_Water_Amp_Phase_FreqMap.ps'];
if exist(s);delete(s);end
figs=figure('visible', 'off');

 plotImage= Vol2Image( Water_rrr );

 imagesc(abs(plotImage),[0 max(abs(Water_rrr(:)))]);%,[ 0, 10*mean(image2plot(:))] )
colormap default;colorbar;title('Water Signal Amplitude')
print(figs, '-append', '-dpsc2', s);

imagesc(angle(plotImage));%,[ 0, 10*mean(image2plot(:))] )
colormap default;colorbar;title('Water Signal Phase')
print(figs, '-append', '-dpsc2', s);

 plotImage= Vol2Image( WaterFreqMap_rrr );
imagesc(plotImage,[-30 30]);%,[ 0, 10*mean(image2plot(:))] )
colormap default;colorbar;title('Freq Map Based on Water Signal')
print(figs, '-append', '-dpsc2', s);


% plotImage= Vol2Image( mrsiReconParams.b0map );
%imagesc(plotImage,[-30 30]);%,[ 0, 10*mean(image2plot(:))] )
%colormap default;colorbar;title('Measured B0 Field map')
%print(figs, '-append', '-dpsc2', s);

close all;


fprintf('Coil profile computation and water reconstruction finished.\n');
end
