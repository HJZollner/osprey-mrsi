function Plot2DSpectra(reconResults)


% ************************************************************************
% Plot Diagonal spectra 
%*************************************************************************

%Reference DataSet

DataOrig_frrr=fft(reconResults.Original_RePhased_Data_trrr,[],1);
Data_frrr=fft(reconResults.Recon_RePhased_Data_trrr,[],1);
DataW_frrr=fftshift(fft(reconResults.Water_trrr,[],1),1);

TimeSize=size(DataOrig_frrr,1);
Final_ppm=reconResults.ppm; %(-4.7+((1:TimeSize)*reconResults.mrProt.samplerate/(TimeSize*reconResults.mrProt.NMRFreq)));
[~,maxppm]=min(abs(-1.5 - reconResults.ppm));
%pts=reconResults.MinPPM_pt:reconResults.MaxPPM_pt;
pts=reconResults.MinPPM_pt:maxppm;

marginx=5;
marginy=1;
xstep=1.1;
jump=1;

s=[reconResults.Log_Dir,filesep,'2DSpectra_plot_Recon', reconResults.NameData, '.ps'];
if exist(s);delete(s);end
figs=figure('units','inches');
for z=1:size(Data_frrr,4)
    ystep=2*max(abs(squeeze(Data_frrr(pts,size(Data_frrr,2)/2,size(Data_frrr,3)/2,z))));% 0.01 Braino %0.1 Invivo
    clf
    hold on;
    for x=(1+marginx):jump:(size(Data_frrr,2)-marginx)
        for y=(1+marginy):jump:(size(Data_frrr,3)-marginy)
            
            plot(xstep*(Final_ppm(pts(end))-Final_ppm(pts(1)))*(y-1)/jump+Final_ppm(pts),abs(squeeze(Data_frrr(pts,x,y,z)))+x/jump*ystep,'g-' )
        end
    end
    
    %axis([Final_ppm(pts(1)) Final_ppm(pts(end))+2*(Final_ppm(pts(end))-Final_ppm(pts(1)))+(Final_ppm(pts(end))-Final_ppm(pts(1)))*size(Data_frr,2)/jump 0 size(Data_frr,2)/jump*ystep])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    pos = get(gcf,'pos');
    set(gcf,'pos',[pos(1) pos(2) 20 20])
    
    print(figs,'-bestfit','-append', '-dpsc2',s );
  
end


s=[reconResults.Log_Dir,filesep,'2DSpectra_plot_Orig', reconResults.NameData, '.ps'];
if exist(s);delete(s);end

ystep=2*max(abs(squeeze(DataOrig_frrr(pts,size(DataOrig_frrr,2)/2,size(DataOrig_frrr,3)/2,size(DataOrig_frrr,4)/2))));% 0.01 Braino %0.1 Invivo

for z=1:size(DataOrig_frrr,4)
   clf
    hold on;
    for x=(1+marginx):jump:(size(DataOrig_frrr,2)-marginx)
        for y=(1+marginy):jump:(size(DataOrig_frrr,3)-marginy)
            plot(xstep*(Final_ppm(pts(end))-Final_ppm(pts(1)))*(y-1)/jump+Final_ppm(pts),abs(squeeze(DataOrig_frrr(pts,x,y,z)))+x/jump*ystep,'g-' )
        end
    end
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    pos = get(gcf,'pos');
    set(gcf,'pos',[pos(1) pos(2) 20 20])
    
    print(figs,'-bestfit','-append', '-dpsc2',s );
   
end

pts=1:size(DataW_frrr,1);

s=[reconResults.Log_Dir,filesep,'2DSpectra_plot_Water', reconResults.NameData, '.ps'];
if exist(s);delete(s);end

ystep=2*max(abs(squeeze(DataW_frrr(pts,size(DataW_frrr,2)/2,size(DataW_frrr,3)/2,size(DataW_frrr,4)/2))));% 0.01 Braino %0.1 Invivo
xstep=1.1;
%{
for z=1:size(DataW_frrr,4)
    clf
    hold on;
   
    for x=(1+marginx):jump:(size(DataW_frrr,2)-marginx)
        for y=(1+marginy):jump:(size(DataW_frrr,3)-marginy)
            plot(xstep*(pts(end)-pts(1))*(y-1)/jump+(pts),abs(squeeze(DataW_frrr(pts,x,y,z)))+x/jump*ystep,'g-' )
        end
    end
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    pos = get(gcf,'pos');
    set(gcf,'pos',[pos(1) pos(2) 20 20])
    
    print(figs,'-bestfit','-append', '-dpsc2',s );

end
%}
close all
end