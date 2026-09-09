function PlotDiagSpectraData(reconResults,Data1_frrr,Data2_frrr,name)


% ************************************************************************
% Plot Diagonal spectra 
%*************************************************************************
if exist(name);delete(name);end


TimeSize=size(Data1_frrr,1);
Final_ppm=reconResults.ppm; %(-4.7+((1:TimeSize)*reconResults.mrProt.samplerate/(TimeSize*reconResults.mrProt.NMRFreq)));

pts=1:size(Data1_frrr,1);%reconResults.MinPPM_pt:reconResults.MaxPPM_pt;

fact=norm(Data1_frrr(:))/norm(Data2_frrr(:));
xstep=0.1*pts(end);
ystep=max(abs(squeeze(Data1_frrr(pts,size(Data1_frrr,2)/2,size(Data1_frrr,2)/2))));% 0.01 Braino %0.1 Invivo
for z=1:size(Data1_frrr,4)
    figs=figure();
    hold on;
    for x=1:size(Data1_frrr,2)
        y=size(Data1_frrr,2)-x+1;
        plot(pts(end)+pts+x*xstep,abs(squeeze(Data1_frrr(pts,x,y,z)))+(size(Data1_frrr,2)-x)*ystep,'g-' )
        plot(pts+x*xstep,fact*abs(squeeze(Data2_frrr(pts,x,y,z)))+(size(Data1_frrr,2)-x)*ystep, 'r-')
        
    end
    
    axis([pts(1) pts(end)+2*(pts(end)-pts(1))+size(Data1_frrr,2)*xstep 0 size(Data1_frrr,2)*ystep])
    title(['Diagonal spectra, z-pos:' num2str(z)])
    print(figs,'-bestfit','-append', '-dpsc2', name);
    
    close all;
end

end