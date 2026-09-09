function VisualizeMasks( BodyWater_ctkkk,HMask, BMask ,SMask,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

s=sprintf('%s_Masks_Images.ps',name_fig);
if exist(s);delete(s);end

BodyWater_crrr=squeeze(ifft(ifft(ifft(mean(BodyWater_ctkkk(:,1:5,:,:,:),2),[],3),[],4),[],5));
mrsiData_rrr = squeeze(sum(abs(BodyWater_crrr).^2,1));

if(numel(mrsiData_rrr(:))~= numel(HMask(:)))
    
    %Make Anatomical masks at the mrsiData size

  
    
    Sdata=size(BodyWater_ctkkk);
    SiMask=size(HMask);
    
    [Xm,Ym,Zm] = meshgrid(linspace(1,Sdata(4),SiMask(2)),linspace(1,Sdata(3),SiMask(1)),linspace(1,Sdata(5),SiMask(3)));
    [Xd,Yd,Zd] = meshgrid(1:Sdata(4),1:Sdata(3),1:Sdata(5));
    clear TempBrainMask TempImMask TempSkMask

    TempBrainMask = interp3(Xm,Ym,Zm,BMask,Xd,Yd,Zd,'linear');
    TempImMask = interp3(Xm,Ym,Zm,HMask,Xd,Yd,Zd,'linear');
    TempSMask = interp3(Xm,Ym,Zm,SMask,Xd,Yd,Zd,'linear');
        
    BMask=round(TempBrainMask/max(TempBrainMask(:)));
    HMask=round(TempImMask/max(TempImMask(:)));
    SMask=round(TempSMask/max(TempSMask(:)));
    
    
    
end

figs=figure('visible', 'off');
 image2plot=Vol2Image(squeeze( permute(mrsiData_rrr,[ 2 1 3])));
subplot(1,2,1),imagesc(image2plot);
axis('off');title(' Water Data unmasked')


print(figs, '-append', '-dpsc2', s);
close all;

figs=figure('visible', 'off');
 image2plot=Vol2Image(squeeze( permute(HMask.*mrsiData_rrr,[ 2 1 3])));
subplot(1,2,1),imagesc(image2plot);
axis('off');title(' Water Data inside mask')
 image2plot=Vol2Image(squeeze( permute(HMask,[ 2 1 3])));
subplot(1,2,2),imagesc(image2plot);
axis('off');title(' Head Mask')

print(figs, '-append', '-dpsc2', s);
close all;

figs=figure('visible', 'off');

 image2plot=Vol2Image(squeeze( permute(BMask.*mrsiData_rrr,[ 2 1 3])));
subplot(1,2,1),imagesc(image2plot);
axis('off');title(' Water Data inside mask')
 image2plot=Vol2Image(squeeze( permute(BMask,[ 2 1 3])));
subplot(1,2,2),imagesc(image2plot);
axis('off');title(' Brain Mask')

print(figs, '-append', '-dpsc2', s);
close all;

figs=figure('visible', 'off');
 image2plot=Vol2Image(squeeze( permute(SMask.*mrsiData_rrr,[ 2 1 3])));
subplot(1,2,1),imagesc(image2plot);
axis('off');title(' Water Data inside mask')
 image2plot=Vol2Image(squeeze( permute(SMask,[ 2 1 3])));
subplot(1,2,2),imagesc(image2plot);
axis('off');title(' Skull Mask')

print(figs, '-append', '-dpsc2', s);
close all;





end

