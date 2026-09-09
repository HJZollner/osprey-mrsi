function SENSE = ComputeRoemersProfiles(mrsiReconParams,Water_ctkkk)
make_figures=0;
reduction = 2^(-8);     % usually there is no need to change this
alpha = 10;%mrsiReconParams.UndersamplingF*mrsiReconParams.mu_tv/reduction; % usually there is no need to change this
maxit = 500;            % use 1000 Iterations for optimal image quality
innerIter = 20;         % usually there is no need to change this
%kmask = mrsiReconParams.kmask ;
%Water_crr=squeeze(sum(ifft(ifft(mrsiReconParams.Water_ctkk(:,1:round(end/10),:,:),[],3),[],4),2));

MSize_data=size(mrsiReconParams.mrsiData_ctkkk);

WSize_data=size(Water_ctkkk);

RowSet=[1:ceil(WSize_data(3)/2), ceil(MSize_data(3)-WSize_data(3)/2+1):MSize_data(3)];
ColSet= [1:ceil(WSize_data(4)/2), ceil(MSize_data(4)-WSize_data(4)/2+1):MSize_data(4)];
SlcSet= [1:ceil(WSize_data(5)/2.0), ceil(MSize_data(5)-WSize_data(5)/2.0+1):MSize_data(5)];

%kmask = 0*mrsiReconParams.kmask;
%kmask(RowSet,ColSet) = mrsiReconParams.kmask ;


%Gaussian Kernel
%sigma=size(mrsiReconParams.Water_ctkk,3)/4;
sigma=size(Water_ctkkk,3)/3; %Optimized with comparison to the Full Res Water Data

[X,Y,Z] = ndgrid(1:MSize_data(3), 1:MSize_data(4),1:MSize_data(5));
xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;zc=floor(MSize_data(5)/2)+1;
exponent = -((X-xc).^2 + (Y-yc).^2+ (Z-zc).^2)./(2*sigma^2);
%Kernel = fftshift(1 / (2 * sqrt(2*pi)).*exp(exponent));
Kernel = fftshift(exp(exponent)); % no need to normalize

% Zero Padding Water and filtering
WaterZeroPad_ctrrr=zeros(WSize_data(1),WSize_data(2), MSize_data(3),MSize_data(4),MSize_data(5));
WaterZeroPad_ctrrr(:,:,RowSet,ColSet,SlcSet)=Water_ctkkk;
Kernel_11kkk=reshape(Kernel,[1,1,size(Kernel)]);
WaterZeroPad_ctrrr  = ifft(ifft(ifft(WaterZeroPad_ctrrr.*Kernel_11kkk,[],3),[],4),[],5);
%WaterZeroPad_ctrr  = ifft(ifft(WaterZeroPad_ctkk,[],3),[],4);
%WaterAmp_crr=squeeze(sum(abs(fft(ifft(ifft(mrsiReconParams.Water_ctkk,[],3),[],4),[],2)),2));%Doesn't work with data with too sharp peak in frequency

clear Kernel_11kkk Kernel exponent

%[Uorig,Sorig,Vorig] = svd(reshape(WaterZeroPad_crrr,WSize_data(1),[]),0);
%US_crrrc=reshape(Uorig*Sorig, WSize_data(1),MSize_data(3),MSize_data(4),MSize_data(5),[]);
%WaterAmp_crr=squeeze(mean(abs(WaterZeroPad_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),2));
%WaterPh_crr=squeeze(angle(sum(WaterZeroPad_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2)));


%%
%{
for k=1:size(mrsiReconParams.Water_ctkk,1)
    %[USreal(:,:,k), funCost] = TGV2_Denoise_MultiMetabVol(im, mrsiReconParams.mu_tv*datanorm);
   SENSE(k,:,:)=exp(1i*squeeze(WaterPh_crr(k,:,:)-WaterPh_crr(mrsiReconParams.mrProt.Main_coil_element,:,:))).*squeeze(WaterAmp_crr(k,:,:)).^(0.5);%./squeeze(sum(WaterAmp_crr,1).^(0.5));
% SENSE(k,:,:)=exp(1i*squeeze(WaterPh_crr(k,:,:))).*squeeze(WaterAmp_crr(k,:,:)).^(0.5);%./squeeze(sum(WaterAmp_crr,1).^(0.5));

end
%}
%SENSE=US_crrrc(:,:,:,:,1);% 1st SVD component represents water at best

SENSE=squeeze(mean(WaterZeroPad_ctrrr(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),2));
clear WaterZeroPad_ctrrr 

Data_name=mrsiReconParams.NameData;
s1=[mrsiReconParams.Log_Dir,filesep,'RoemersProfiles_',Data_name,'.ps'];
s2=[mrsiReconParams.Log_Dir,filesep,'RoemersOrigSignal_',Data_name,'.ps'];
if exist(s1);delete(s1);end
if exist(s2);delete(s2);end
figs=figure('visible', 'off');

SENSESQ = sum(abs(SENSE).^2,1);
for k=1:size(Water_ctkkk,1)
    SENSE(k,:,:,:)=SENSE(k,:,:,:)./SENSESQ;
    if make_figures
        subplot(2,3,1),imagesc(abs(squeeze(SENSE(k,:,:,round(end/2)))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,2),imagesc(abs(squeeze(SENSE(k,:,round(end/2),:))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,3),imagesc(abs(squeeze(SENSE(k,round(end/2),:,:))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,4),imagesc(angle(squeeze(SENSE(k,:,:,round(end/2)))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,5),imagesc(angle(squeeze(SENSE(k,:,round(end/2),:))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,6),imagesc(angle(squeeze(SENSE(k,round(end/2),:,:))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        print(figs, '-append', '-dpsc2', s1);
        
        subplot(2,3,1),imagesc(abs(squeeze(US_crrrc(k,:,:,round(end/2),1))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,2),imagesc(abs(squeeze(US_crrrc(k,:,round(end/2),:,1))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,3),imagesc(abs(squeeze(US_crrrc(k,round(end/2),:,:,1))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,4),imagesc(angle(squeeze(US_crrrc(k,:,:,round(end/2),1))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,5),imagesc(angle(squeeze(US_crrrc(k,:,round(end/2),:,1))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        subplot(2,3,6),imagesc(angle(squeeze(US_crrrc(k,round(end/2),:,:,1))));%,[ 0, 10*mean(image2plot(:))] )
        colormap default;
        print(figs, '-append', '-dpsc2', s2);
    end
end
SENSE(find(isnan(SENSE)))=0;
clear Water_crrr
close all;
end
