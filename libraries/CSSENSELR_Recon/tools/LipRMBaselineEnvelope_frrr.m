function  Corrected_Data_frrr    = LipRMBaselineEnvelope_frrr( mrsiData_frrr,TE,Mask,mrsiReconParams,NameData)

NbT=size(mrsiData_frrr,1);
MaxPPM_BLCorr=-1.5;
SampleRate = mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;
t=(1:NbT)/SampleRate;


FreqOI=mrsiReconParams.MinPPM_pt:mrsiReconParams.MaxPPM_pt;
%FreqOI=ones(NbT,1);

if nargin>=5
    if NbT>(mrsiReconParams.MaxPPM_pt-mrsiReconParams.MinPPM_pt+1)
    	FreqPlot=mrsiReconParams.MinPPM_pt:mrsiReconParams.MaxPPM_pt;
    	else
    	FreqPlot=1:NbT;
    end
    
    s1=[mrsiReconParams.Log_Dir,filesep,NameData,'_SVD_SpecComp_BeforeBLRM.ps'];
    s2=[mrsiReconParams.Log_Dir,filesep,NameData,'_SVD_SpatComp_BeforeBLRM.ps'];
    delete(s1);delete(s2);
    mrsiData_rrrf=permute(mrsiData_frrr,[2,3,4,1]);

    SizeD=size(mrsiData_rrrf);
    [Uorig,Sorig,Vorig] = svd(reshape(mrsiData_rrrf(:,:,:,FreqPlot).*mrsiReconParams.BrainMask,[],numel(FreqPlot)),0);
    U_rrrc=reshape(Uorig(:,1:10), SizeD(1),SizeD(2),SizeD(3),[]);
    figs=figure('visible', 'off');
    for comp=1:10
        
        plot(FreqPlot,real(Vorig(:,comp)),...
            FreqPlot,imag(Vorig(:,comp)),...
            FreqPlot,abs(Vorig(:,comp)));
        title(['Spectral comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s1);
        
        plotImage= Vol2Image(abs(U_rrrc(:,:,:,comp)) );
        imagesc(plotImage);
        colormap default;colorbar;
        title(['Spatial comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s2);
        
    end
    clear MRSIDataLR_rrrf
end
close all

T=1:round(25E-3*SampleRate);% only the 25 last ms of the FID
Apod_t111=[ones(1,NbT-numel(T)) exp(-(T.^2)*(4/numel(T)).^2)];
Apod_t111=reshape(Apod_t111,[NbT,1,1,1]);
FOPhase=2*pi/NbT*TE*SampleRate;

Mask_tr=repmat(Mask(:)',[NbT,1]);

mrsiData_fr = mrsiData_frrr(Mask_tr>0);
mrsiData_fr = reshape(mrsiData_fr,NbT,[]);

minNbPt=round(NbT*30*(mrsiReconParams.FieldStrength/7.0)/mrsiReconParams.mrProt.samplerate);
maxNbPt=round(NbT*60*(mrsiReconParams.FieldStrength/7.0)/mrsiReconParams.mrProt.samplerate);

%minNbPt=round(NbT*100*(mrsiReconParams.FieldStrength/7.0)/mrsiReconParams.mrProt.samplerate);
%maxNbPt=round(NbT*200*(mrsiReconParams.FieldStrength/7.0)/mrsiReconParams.mrProt.samplerate);

MinPPM_pt=mrsiReconParams.MinPPM_pt;
MaxPPM_pt=mrsiReconParams.MaxPPM_pt;

%{
[~,WinStart_pt]=min(abs((mrsiReconParams.MinPPM+0.5)+mrsiReconParams.PPMshift  - mrsiReconParams.ppm));
[~,WinEnd_pt]=min(abs((MaxPPM_BLCorr)+mrsiReconParams.PPMshift  - mrsiReconParams.ppm));


FrSt=1:round(WinStart_pt-MinPPM_pt);
FrEnd=1:round(MaxPPM_pt-WinEnd_pt);
Apod_f=[ flip(exp(-(FrSt.^2)*(2/numel(FrSt)).^2)) ones(1,size(Short_R_BL_fr,1)-numel(FrEnd)-numel(FrSt)) exp(-(FrEnd.^2)*(2/numel(FrEnd)).^2) zeros(1,NbT - size(Short_R_BL_fr,1))]';
%}

Apod_f=[exp(-((mrsiReconParams.ppm-(-2.0))/0.5).^2) ]';

Short_R_BL_fr=0*mrsiData_fr(MinPPM_pt:MaxPPM_pt,:);
HalfRange= round(0.33*(MaxPPM_pt-MinPPM_pt));
parfor a=1:size(mrsiData_fr,2)
    if norm(mrsiData_fr(:,a))>0
  	[ mrsiData_fr(:,a) ,ZOPhase_r(a)] = Rephase_poly(mrsiData_fr(:,a).*exp(-1j*((1:NbT)'*FOPhase))  );
    Spectrum=mrsiData_fr(:,a);
    Spectrum=Spectrum(MinPPM_pt:MaxPPM_pt);
    MVal=max(real(Spectrum(:)));
    PaddSpectrum=[randn(HalfRange,1)*1E-2*MVal; Spectrum ; randn(HalfRange,1)*1E-2*MVal];
    BL=zeros(numel(PaddSpectrum),numel(minNbPt:maxNbPt));
    for NbPtBL=minNbPt:maxNbPt
        [~,BL(:,(NbPtBL-minNbPt+1))] = envelope(real(double(PaddSpectrum)),NbPtBL,'peak'); 
        %[~,BL(:,(NbPtBL-minNbPt+1))] = envelope(real(double(PaddSpectrum))); 
        %BL(:,(NbPtBL-minNbPt+1))=smooth(BL(:,(NbPtBL-minNbPt+1)),NbPtBL,'lowess');
    end

    Short_R_BL_fr(:,a)=mean(BL((HalfRange+1):(end-HalfRange),:),2);

    else
        ZOPhase_r(a) =0;
    end
end

mrsiData_tr = ifft(mrsiData_fr,[],1);
BL_fr=0*mrsiData_fr;
BL_fr(MinPPM_pt:MaxPPM_pt,:)=Short_R_BL_fr+1j*Compute_ImagPart_From_KK(Short_R_BL_fr);;
%BL_fr(MinPPM_pt:MaxPPM_pt,:)=Short_R_BL_fr;
BL_fr=BL_fr.*Apod_f;
BL_fr=transpose(transpose(BL_fr)*(eye(NbT)-mrsiReconParams.LipidRM_ff));
Corrected_Data_tr=ifft(mrsiData_fr-BL_fr,[],1);


BL_frrr = mrsiData_frrr;
BL_frrr = reshape(BL_frrr,NbT,[]);
BL_fr=fft(mrsiData_tr - Corrected_Data_tr,[],1);
BL_frrr(Mask_tr>0)=BL_fr.*exp(1j*(ZOPhase_r + (1:NbT)'*FOPhase));
BL_frrr= reshape(fft(ifft(BL_frrr,[],1).*Apod_t111,[],1),size(mrsiData_frrr));

Corrected_Data_frrr = mrsiData_frrr - BL_frrr;
 
if nargin>=5
    MRSIDataLR_rrrf=permute(Corrected_Data_frrr,[2,3,4,1]);
    MRSIDataOrig_rrrf=permute(mrsiData_frrr,[2,3,4,1]);

    
    s=[mrsiReconParams.Log_Dir,filesep,NameData,'_Diagonal_Spectra_BeforeAndAfterBLRM.ps'];
    delete(s);  
    figs=figure('visible', 'off');   
	for c=1:4:size(MRSIDataLR_rrrf,3);
		for a=1:4:min(size(MRSIDataLR_rrrf,1),size(MRSIDataLR_rrrf,2));if(mrsiReconParams.BrainMask(a,a,c)==1)
            
                subplot(2,2,1);title('Real part')
		
                sp1=squeeze(MRSIDataLR_rrrf(a,a,c,FreqPlot));
		sp2=squeeze(MRSIDataOrig_rrrf(a,a,c,FreqPlot));

    		 plot(FreqPlot,real(sp1),...
           	 FreqPlot,real(sp2),...
            	FreqPlot,real(sp1-sp2));

		subplot(2,2,2);title('Imag part')
    		 plot(FreqPlot,imag(sp1),...
           	 FreqPlot,imag(sp2),...
            	FreqPlot,imag(sp1-sp2));

		subplot(2,2,3);title('Magnitude')
    		 plot(FreqPlot,abs(sp1),...
           	 FreqPlot,abs(sp2),...
            	FreqPlot,abs(sp1-sp2));
                title(['x=',num2str(a) , 'y=',num2str(a) ,'z=',num2str(c)]);
                %legend({'Corrected','Original','Difference'},'Location','eastoutside')
 		subplot(2,2,4);
            plot(FreqPlot,ones(size(FreqPlot)),...
           	 FreqPlot,ones(size(FreqPlot)),...
            	FreqPlot,ones(size(FreqPlot)));
               legend({'Corrected','Original','Difference'})

    		print(figs, '-append', '-dpsc2', s);
    		clf
    		
		end;end
	end
    clear MRSIDataOrig_rrrf
    close all;
    
    s1=[mrsiReconParams.Log_Dir,filesep,NameData,'_SVD_SpecComp_AfterBLRM.ps'];
    s2=[mrsiReconParams.Log_Dir,filesep,NameData,'_SVD_SpatComp_AfterBLRM.ps'];
    delete(s1);delete(s2);
	
    [Uorig,Sorig,Vorig] = svd(reshape(MRSIDataLR_rrrf(:,:,:,FreqPlot).*mrsiReconParams.BrainMask,[],numel(FreqPlot)),0);
    clear MRSIDataOrig_rrrf
    U_rrrc=reshape(Uorig(:,1:10), SizeD(1),SizeD(2),SizeD(3),[]);
    figs=figure('visible', 'off');
    for comp=1:10
        
        plot(FreqPlot,real(Vorig(:,comp)),...
            FreqPlot,imag(Vorig(:,comp)),...
            FreqPlot,abs(Vorig(:,comp)));
        title(['Spectral comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s1);
        
        plotImage= Vol2Image(abs(U_rrrc(:,:,:,comp)) );
        imagesc(plotImage);
        colormap default;colorbar;
        title(['Spatial comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s2);
        
    end
end

end

