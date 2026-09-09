function VisualizeTGV( Vol_After ,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here
warning('off','MATLAB:DELETE:FileNotFound');
nDims=ndims(Vol_After);
SizeData=size(Vol_After);
s=sprintf('%s_SpatialComp.ps',name_fig);
delete(s); 
%ImSiC=2048;ImSiR=2048;
for k=1:SizeData(end);

        figs(k) = figure('visible', 'off'); %,'Position', [100, 200, ImSiC,ImSiR]); 
      
        
        image2plot = Vol2Image( permute(Vol_After(1:(end),1:(end),:,k),[ 2 1 3 4] ) );;

        lim = quantile(abs(image2plot(:)),0.98);
       subplot(1,2,1);  imagesc(abs(image2plot),[0 lim]), title('Recon Combination Amp.'),daspect([1 1 1]);  
        axis('off');
        colormap default 
        
        subplot(1,2,2);  imagesc(angle(image2plot)), title('Recon Combination Phase.'),daspect([1 1 1]);  
        axis('off');
        colormap default 
        
        orient(figs(k),'landscape');
        print(figs(k),'-bestfit', '-append', '-dpsc2', s); 
        
end;
close all;




end

