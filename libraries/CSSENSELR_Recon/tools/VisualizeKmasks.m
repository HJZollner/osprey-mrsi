function VisualizeKmasks( mrsiData_ctkkk,kmask_data, HeadWater_ctkkk ,kmask_Head,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

s=sprintf('%s_Kmasks_Images.ps',name_fig);
if exist(s);delete(s);end
Size_data=size(mrsiData_ctkkk);

kmask_data=fftshift(kmask_data);
kmask_Head=fftshift(kmask_Head);

mrsiData_kkk = log(squeeze(sum(sum(abs(fft(mrsiData_ctkkk,[],2)),1),2)));
HWater_kkk = log(squeeze(sum(sum(abs(fft(HeadWater_ctkkk,[],2)),1),2)));
mrsiData_kkk=fftshift(mrsiData_kkk);
HWater_kkk=fftshift(HWater_kkk);

figs=figure('visible', 'off');
subplot(2,2,1),imagesc(Vol2Image(mrsiData_kkk));
axis('off');
colormap default ;title(' Log MRSI Data k-space Combined Amplitude')
subplot(2,2,3),imagesc(Vol2Image(kmask_data));
axis('off');
colormap default; title(' MRSI Data k-mask')
subplot(2,2,2),imagesc(Vol2Image(HWater_kkk));
axis('off');
colormap default; title(' Log Head Water k-space Combined Amplitude')
subplot(2,2,4),imagesc(Vol2Image(kmask_Head));
axis('off');
colormap default ; title(' Head Water k-mask')
print(figs, '-append', '-dpsc2', s);
close all;

figs=figure('visible', 'off');
subplot(2,2,1),imagesc(mrsiData_kkk(:,:,round(end/2)));
axis('off');
colormap default ;title('Axial Log MRSI Data k-space Combined Amplitude')
subplot(2,2,3),imagesc(kmask_data(:,:,round(end/2)));
axis('off');
colormap default; title('Axial MRSI Data k-mask')
subplot(2,2,2),imagesc(HWater_kkk(:,:,round(end/2)));
axis('off');
colormap default; title('Axial Log Head Water k-space Combined Amplitude')
subplot(2,2,4),imagesc(kmask_Head(:,:,round(end/2)));
axis('off');
colormap default ; title('Axial Head Water k-mask')
print(figs, '-append', '-dpsc2', s);
close all;

figs=figure('visible', 'off');
subplot(2,2,1),imagesc(squeeze(mrsiData_kkk(:,round(end/2),:)));
axis('off');
colormap default ;title('Coro Log MRSI Data k-space Combined Amplitude')
subplot(2,2,3),imagesc(squeeze(kmask_data(:,round(end/2),:)));
axis('off');
colormap default; title('Coro MRSI Data k-mask')
subplot(2,2,2),imagesc(squeeze(HWater_kkk(:,round(end/2),:)));
axis('off');
colormap default; title('Coro Log Head Water k-space Combined Amplitude')
subplot(2,2,4),imagesc(squeeze(kmask_Head(:,round(end/2),:)));
axis('off');
colormap default ; title('Coro Head Water k-mask')
print(figs, '-append', '-dpsc2', s);
close all;

figs=figure('visible', 'off');
subplot(2,2,1),imagesc(squeeze(mrsiData_kkk(round(end/2),:,:)));
axis('off');
colormap default ;title('Sag Log MRSI Data k-space Combined Amplitude')
subplot(2,2,3),imagesc(squeeze(kmask_data(round(end/2),:,:)));
axis('off');
colormap default; title('Sag MRSI Data k-mask')
subplot(2,2,2),imagesc(squeeze(HWater_kkk(round(end/2),:,:)));
axis('off');
colormap default; title('Sag Head Log Water k-space Combined Amplitude')
subplot(2,2,4),imagesc(squeeze(kmask_Head(round(end/2),:,:)));
axis('off');
colormap default ; title('Sag Head Water k-mask')
print(figs, '-append', '-dpsc2', s);
close all;




end

