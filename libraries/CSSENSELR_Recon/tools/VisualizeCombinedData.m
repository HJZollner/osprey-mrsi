function VisualizeCombinedData( mrsiData, HeadWater_rrr ,mrsiReconParams,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

s=sprintf('%s_Amp_CoilCombined_3DRawData_Images.ps',name_fig);
if exist(s);delete(s);end
s2=sprintf('%s_Phase_CoilCombined_3DRawData_Images.ps',name_fig);
if exist(s2);delete(s2);end
Size_data=size(mrsiData);


Sdata=size(HeadWater_rrr);
SMask=size(squeeze(mrsiReconParams.SENSE(1,:,:,:)));
    
%[Xm,Ym,Zm] = meshgrid(linspace(1,Sdata(4),SMask(2)),linspace(1,Sdata(3),SMask(1)),linspace(1,Sdata(5),SMask(3)));
%[Xd,Yd,Zd] = meshgrid(1:Sdata(4),1:Sdata(3),1:Sdata(5));
%for c=1:size(mrsiReconParams.SENSE,1)
%	WSENSE(c,:,:,:) = interp3(Xm,Ym,Zm,squeeze(mrsiReconParams.SENSE(c,:,:,:)),Xd,Yd,Zd,'spline');
%end

%WImMask=interp3(Xm,Ym,Zm,mrsiReconParams.ImMask,Xd,Yd,Zd,'nearest');

SENSE=mrsiReconParams.SENSE./sum(abs(mrsiReconParams.SENSE).^2,1);
%WSENSE=WSENSE./sum(abs(WSENSE).^2,1);


%HWater_ctrrr = ifft(ifft(ifft(HeadWater,[],3),[],4),[],5);

%WaterAmp_crrr=squeeze(mean(abs(fft(HWater_ctrrr,[],2)),2));
%WaterPh_crrr=squeeze(angle(sum(HWater_ctrrr(:,1:5,:,:,:),2)));
%Head_crrr=WaterAmp_crrr.*exp(1j*WaterPh_crrr);
%Head_crrr=squeeze(mean(HWater_ctrrr(:,1:5,:,:,:),2));

%CombHead_rrr=squeeze(sum(Head_crrr.*conj(WSENSE),1)).*WImMask;
%clear HWater_ctrrr  WaterAmp_crrr WaterPh_crrr
CombHead_rrr=HeadWater_rrr;

MRSIData_ctrrr =ifft(ifft(ifft(mrsiData,[],3),[],4),[],5);

%WaterAmp_crrr=squeeze(mean(abs(fft(MRSIData_ctrrr,[],2)),2));
%WaterPh_crrr=squeeze(angle(sum(MRSIData_ctrrr(:,1:5,:,:,:),2)));
%Data_crrr=WaterAmp_crrr.*exp(1j*WaterPh_crrr);
Data_crrr=squeeze(mean(MRSIData_ctrrr(:,1:5,:,:,:),2));
CombData_rrr=squeeze(sum(Data_crrr.*conj(SENSE),1));

clear MRSIData_ctrrr WaterAmp_crrr  WaterPh_crrr


   figs=figure('visible', 'off'); 
         
       image2plot=Vol2Image(abs(CombData_rrr));
       subplot(1,2,1);imagesc(abs(image2plot));
       title('MRSI Data Amp.'); axis('off');colormap default ;
       
        image2plot=Vol2Image(abs(CombHead_rrr));
       subplot(1,2,2);imagesc(abs(image2plot)),
       title('Head Water Data Amp.'); axis('off');colormap default ;
        
        print(figs, '-append', '-dpsc2', s); 

 	image2plot=Vol2Image(abs(CombData_rrr).*mrsiReconParams.ImMask);
       subplot(1,2,1);imagesc(abs(image2plot));
       title('MRSI Data Amp. Head Mask'); axis('off');colormap default ;
       
        image2plot=Vol2Image(abs(CombHead_rrr).*mrsiReconParams.ImMask);
       subplot(1,2,2);imagesc(abs(image2plot)),
       title('Head Water Data Amp. Head Mask'); axis('off');colormap default ;
        
        print(figs, '-append', '-dpsc2', s); 

 	image2plot=Vol2Image(abs(CombData_rrr).*mrsiReconParams.BrainMask);
       subplot(1,2,1);imagesc(abs(image2plot));
       title('MRSI Data Amp. Brain Mask'); axis('off');colormap default ;
       
        image2plot=Vol2Image(abs(CombHead_rrr).*mrsiReconParams.BrainMask);
       subplot(1,2,2);imagesc(abs(image2plot)),
       title('Head Water Data Amp. Brain Mask'); axis('off');colormap default ;
        
        print(figs, '-append', '-dpsc2', s); 

	image2plot=Vol2Image(abs(CombData_rrr).*mrsiReconParams.SkMask);
       subplot(1,2,1);imagesc(abs(image2plot));
       title('MRSI Data Amp. Skull Mask'); axis('off');colormap default ;
       
        image2plot=Vol2Image(abs(CombHead_rrr).*mrsiReconParams.SkMask);
       subplot(1,2,2);imagesc(abs(image2plot)),
       title('Head Water Data Amp. Skull Mask'); axis('off');colormap default ;
        
        print(figs, '-append', '-dpsc2', s); 

close all;


   figs=figure('visible', 'off'); 
   
        image2plot=Vol2Image(angle(CombData_rrr));
       subplot(1,2,1);imagesc(angle(image2plot));
       title('MRSI Data Phase'); axis('off');colormap default ;
       
        image2plot=Vol2Image(angle(CombHead_rrr));
       subplot(1,2,2);imagesc(angle(image2plot)),
       title('Head Water Data Phase'); axis('off');colormap default ;
        
        
        print(figs, '-append', '-dpsc2', s2); 

close all;


end

