function VisualizeData( mrsiData_ctkkk, HeadWater_ctkkk ,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

s=sprintf('%s_Amp_3DRawData_Images.ps',name_fig);
if exist(s);delete(s);end
s2=sprintf('%s_Phase_3DRawData_Images.ps',name_fig);
if exist(s2);delete(s2);end
Size_data=size(mrsiData_ctkkk);

Head_crrr =squeeze(ifft(ifft(ifft(mean(HeadWater_ctkkk(:,1:5,:,:,:),2),[],3),[],4),[],5));


clear HWater_ctrrr  WaterAmp_crrr WaterPh_crrr
Data_crrr =squeeze(ifft(ifft(ifft(mean(mrsiData_ctkkk(:,1:5,:,:,:),2),[],3),[],4),[],5));


for C=1:Size_data(1);
   figs(C)=figure('visible', 'off'); 
         
       image2plot=Vol2Image(squeeze( permute(Data_crrr(C,1:end,1:end,:),[ 1 3 2 4])));
       subplot(1,2,1);imagesc(abs(image2plot));
       title('MRSI Data Amp.'); axis('off');colormap default ;
       
        image2plot=Vol2Image(squeeze(permute(Head_crrr(C,1:end,1:end,:),[ 1 3 2 4])));
       subplot(1,2,2);imagesc(abs(image2plot)),
       title('Head Water Data Amp.'); axis('off');colormap default ;
        
        print(figs(C), '-append', '-dpsc2', s); 
    end;
close all;

for C=1:Size_data(1);
   figs(C)=figure('visible', 'off'); 
   
        image2plot=Vol2Image(squeeze(permute( Data_crrr(C,1:end,1:end,:),[ 1 3 2 4])));
       subplot(1,2,1);imagesc(angle(image2plot));
       title('MRSI Data Amp.'); axis('off');colormap default ;
       
        image2plot=Vol2Image(squeeze(permute(Head_crrr(C,1:end,1:end,:),[ 1 3 2 4])));
       subplot(1,2,2);imagesc(angle(image2plot)),
       title('Head Water Data Amp.'); axis('off');colormap default ;
        
        
        print(figs(C), '-append', '-dpsc2', s2); 
    end;
close all;


end

