function PlotDiagSpectra(reconResults)


% ************************************************************************
% Plot Diagonal spectra 
%*************************************************************************

%Reference DataSet
if isfield(reconResults,'mrsiData_wLip_trrr')
    DatawL_frrr=fft(reconResults.mrsiData_wLip_trrr,[],1);
elseif isfield(reconResults,'Original_RePhased_Data_trrr')
    DatawL_frrr=fft(reconResults.Original_RePhased_Data_trrr,[],1);
else
    error('No reference dataset found in reconResults!')
end

%Resulting DataSet
if isfield(reconResults,'mrsiDataLipRem_trrr')
  Data_frrr=fft(reconResults.mrsiDataLipRem_trrr,[],1);
elseif isfield(reconResults,'Recon_RePhased_Data_trrr')
    Data_frrr=fft(reconResults.Recon_RePhased_Data_trrr,[],1);
 else
    error('No resulting dataset found in reconResults!')
end   

%2nd Resulting DataSet    
if isfield(reconResults,'mrsiDataLipRem2_trrr')
    Data2_frrr=fft(reconResults.mrsiDataLipRem2_trrr,[],1);
end

TimeSize=size(DatawL_frrr,1);
Final_ppm=reconResults.ppm; %(-4.7+((1:TimeSize)*reconResults.mrProt.samplerate/(TimeSize*reconResults.mrProt.NMRFreq)));

pts=reconResults.MinPPM_pt:reconResults.MaxPPM_pt;

xstep=0.3;
ystep=max(abs(squeeze(Data_frrr(pts,size(Data_frrr,2)/2,size(Data_frrr,2)/2))));% 0.01 Braino %0.1 Invivo
figs=figure(); 
hold on;
for x=1:size(Data_frrr,2)
    y=size(Data_frrr,2)-x+1;
    plot(Final_ppm(pts(end))-Final_ppm(pts(1))+Final_ppm(pts)+x*xstep,abs(squeeze(Data_frrr(pts,x,y)))+(size(Data_frrr,2)-x)*ystep,'g-' )
    plot(Final_ppm(pts)+x*xstep,abs(squeeze(DatawL_frrr(pts,x,y)))+(size(Data_frrr,2)-x)*ystep, 'r-')
    if isfield(reconResults,'mrsiDataLipRem2_trrr')
        plot(2*Final_ppm(pts(end))-2*Final_ppm(pts(1))+Final_ppm(pts)+x*xstep,abs(squeeze(Data2_frrr(pts,x,y)))+(size(Data_frrr,2)-x)*ystep, 'b-')
    end
end 

axis([Final_ppm(pts(1)) Final_ppm(pts(end))+2*(Final_ppm(pts(end))-Final_ppm(pts(1)))+size(Data_frr,2)*xstep 0 size(Data_frr,2)*ystep])

 print(figs, '-dpsc2', [reconResults.Log_Dir,filesep,'DiagonalSpectra_plot_', reconResults.NameData, '.ps']);
 close all;
 
end