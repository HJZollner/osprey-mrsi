clear ;
close all;

%filename{1}='/home/antoine/Brain/MRS/Prisma_data/Data-2016-02-16-AK/MRSI_CFLr_02-16.mat';
%filename{2}='/home/antoine/Brain/MRS/Prisma_data/Data-2016-01-19-AK/Original_TWIX/MRSI_CFLr_01-19.mat';
%filename{3}='/home/antoine/Brain/MRS/Prisma_data/Data-2016-02-22-AK/MRSI_CFLr_02-22.mat';
%filename{4}='/home/antoine/Brain/MRS/Prisma_data/Data-2016-02-19-AK/MRSI_CFLr_02-19.mat';

filename{1}='/home/antoine/Brain/MRS/Prisma_data/Data-2016-02-16-AK/meas_MID00025_FID21908_CSI_ADI_SE_64x64_8ovs_WWS50hz_228x228x10mm_TE30_TR1600.dat';
filename{2}='/home/antoine/Brain/MRS/Prisma_data/Data-2016-01-19-AK/Original_TWIX/meas_MID00482_FID12398_CSI_ADI_SE_2D_64x64_8OVS_220IR_WWS.dat';
filename{3}='/home/antoine/Brain/MRS/Prisma_data/Data-2016-02-22-AK/meas_MID00172_FID23960_CSI_ADI_SE_2D_64x64_8OVS_TE30_TR1600_WWS_220IR.dat';
filename{4}='/home/antoine/Brain/MRS/Prisma_data/Data-2016-02-19-AK/meas_MID00340_FID23791_CSI_ADI_SE_2D_64x64_8OVS_TE30_TR1600_WWS_220IR.dat';


index=1;
for fi=1:numel(filename);
    fi
    %load(filename{fi});
    
     [twix]=reading_twix_small(filename{fi});
    SE=squeeze(sum(sum(abs(twix.raw_tckk).^2,1),2));
    kmask=squeeze((sum(sum(abs(twix.raw_tckk),1),2)>0));
    %imagesc(SE);
    %imagesc(squeeze(abs(sum(fftshift(fftshift(mrsiReconParams.mrsiData,2),3),1))));
    %KData_center=fftshift(fftshift(mrsiReconParams.mrsiData,2),3);
  
    %TimeS_area(:,:,fi)=squeeze(sum(abs(KData_center(1:100,:,:)),1))-squeeze(sum(abs(KData_center((end-100):end,:,:)),1));
    %Spectrum_E(:,:,fi)=squeeze(sum(abs(fft(KData_center,[],1)).^2,1));
    Spectrum_E(:,:,fi)=SE;
    Size_grid=size(Spectrum_E);
 starting_radius=1;;
 
    for a=1:Size_grid(1);
        for b=1:Size_grid(2);
            radius(a,b)=sqrt((a-Size_grid(1)*0.5-1)^2+ (b-Size_grid(2)*0.5-1)^2);
            if(radius(a,b)<=(Size_grid(1)*0.5-1) & radius(a,b)>=starting_radius)
              Data_SpE(index)=Spectrum_E(a,b,fi);
             % Data_SpE(index)=TimeS_area(a,b,fi);
              Data_Radius(index)=radius(a,b);
            index=index+1;
          end
        end
    end
  
    
end
%%
figure();
subplot(2,2,1);
scatter((Data_Radius*64/228)',Data_SpE', 'SizeData', 5)
set(gca,'xscale','log')
set(gca,'Yscale','log')
xlabel('K radius [1/mm]');
ylabel('A.U.')
hold on;
f = fit(log(Data_Radius*64/228)',log(Data_SpE)','poly1');

expo=f.p1
CI=confint(f);
expo_int=(CI(2,1)-CI(1,1))*0.5;

nump=(1*64/228):0.1:(32*64/228);
fitting=(exp(f.p2)*nump.^f.p1);
plot(nump,fitting,'r')
legend('Spectral Density',sprintf('Fitted exponent = %g %c %g ',f.p1,char(177),expo_int))
hold off;

Prob_Map=(radius.^(expo/2.0));
%Prob_Map=(radius.^(mean(expo)));
Fil_Fact=0.25;
ploti=2;
subplot(2,2,ploti);
ploti=ploti+1;
imagesc(kmask);
title(sprintf('K-space sampling with filling factor = %g',1 ));
xlabel('Kx');
ylabel('Ky')

for Fil_Fact=[0.5 0.25];
    Fil_Fact
   % New_SamplingMask=fftshift(fftshift(kmask,1),2);
    New_SamplingMask=kmask;
    Tot_AcqP=sum(kmask(:));
    fil=1;
    while(fil>Fil_Fact)
        a=round(rand()*(Size_grid(1)-1))+1;
        b=round(rand()*(Size_grid(2)-1))+1;
        if(Prob_Map(a,b)*8<rand()) ;
            New_SamplingMask(a,b)=0;  
        end
        fil=sum(New_SamplingMask(:))/Tot_AcqP;
    end
    subplot(2,2,ploti);
    ploti=ploti+1;
    imagesc(New_SamplingMask)
    title(sprintf('K-space sampling with filling factor = %g',Fil_Fact ));
    xlabel('Kx');
ylabel('Ky')

end
map= [ 1 ,1,1; 0,0,1];
colormap(map);
%figure()
%imagesc(New_SamplingMask)
%imagesc(Prob_Map)
%%
for rad_i=1:321
    rad=8+rad_i/10.0;
    Data(rad_i)=0;
    NB_inside=0;
    for a=1:Size_grid(1);
        for b=1:Size_grid(2);
            if sqrt((a-Size_grid(1)*0.5-1)^2+ (b-Size_grid(2)*0.5-1)^2)<rad;
                Data(rad_i)=Data(rad_i)+New_SamplingMask(a,b);
                NB_inside=NB_inside+1;
            end
        end
    end
    Data(rad_i)=Data(rad_i)/NB_inside;
end
