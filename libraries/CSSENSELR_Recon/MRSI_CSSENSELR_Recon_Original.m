function MRSI_CSSENSELR_Recon(MRSI_Filename, Water_Filename, Data_name,FieldStrength, mu_tv, LipidRemStrength, WaterSuppComp, Comp,  NbCoilCompression,UndersamplingF )

% (c) 1.9.2023
% Antoine Klauser (antoine.klauser@gmail.com)
% If you consider this code to be useful for your research, please cite:
%  Klauser A, Klauser P, Grouiller F, Courvoisier S, Lazeyras F. Whole‐brain high‐resolution metabolite mapping with 3D compressed‐sensing SENSE low‐rank 1 H FID‐MRSI.
%  NMR Biomed. 2022;35(1). doi:10.1002/nbm.4615
%
% The TGV reconstruction was build on the TGV recon for MRI provided by Florian Knoll
% Don't forget to mention his work if your work focus on the regularization part :
% Knoll, F.; Bredies, K.; Pock, T.; Stollberger, R.: Second Order Total
% Generalized Variation (TGV) for MRI: Magnetic Resonance in Medicine, 
% to appear (2010)
%
% Reconstruction function inputs:
% MRSI_Filename : the MRSI raw data file name
% Water_Filename : the Water raw data file name
% Data_name : the name given to  the resulting dataset
% FieldStrength : the fieldstrength (3/7 Tesla),
% mu_tv : the TGV weight (1E-4 - 1E-2)
% LipidRemStrength : Lipid suppression strength (0-2, 1 is standard),
% WaterSuppComp : the number of components for water suppression (8-32, 0 for no water suppression), 
% Comp : the number of components for the LowRank decomposition (8-40, depending on the application),
% NbCoilCompression : Coil Compression Nb (99 for none, 12 is safe), 
% UndersamplingF : the undersampling ratio (1 for none)."


[SCRIPT_DIR, ~, ~]=fileparts(mfilename('fullpath'));
addpath(genpath (SCRIPT_DIR));
CURRENT_DIR=pwd;
warning('off','MATLAB:DELETE:FileNotFound')

mrsiReconParams=struct('NameData',Data_name,'mu_tv',mu_tv,'modelOrder',Comp,'MRSI_Filename',MRSI_Filename,'Water_Filename',Water_Filename,'NbCoilCompression',NbCoilCompression);

mrsiReconParams.Log_Dir = ['Log_Files_', Data_name];
mrsiReconParams.Results_Dir = ['Results_',Data_name];
mrsiReconParams.UndersamplingF = UndersamplingF;

mrsiReconParams.AcqDelay = 0.9E-3; %TE ms
mrsiReconParams.SpinEcho = 0; %If Spin Echo, mrsiReconParams.AcqDelay disregarded

mrsiReconParams.FieldStrength=FieldStrength;


mrsiReconParams.FreqWindow4Recon=1; %yes or no : perform the recon on a small frequency window only (improve the performance)?
if mrsiReconParams.FreqWindow4Recon
    if mrsiReconParams.FieldStrength==7
        mrsiReconParams.MinPPM   =  -4.1 ;% start of the frequency window in ppm (downfield) with negative sign
    else
        mrsiReconParams.MinPPM   = -4.3 ;% start of the frequency window in ppm (downfield) with negative sign
    end
    mrsiReconParams.MaxPPM   = -1.0 ; % end of the frequency window in ppm (upfield) with negative sign
else
    mrsiReconParams.MinPPM   =  -999 ;  % No windowing with the following values
    mrsiReconParams.MaxPPM   = 999 ;
end

mrsiReconParams.LimitAcqBW=1; %yes or no : reduce the acquistion bandwidth by interpolation of the FID points? (improve recon performance if bandwith is large)
if mrsiReconParams.LimitAcqBW
    if mrsiReconParams.FieldStrength==7
        mrsiReconParams.MaxSweepwidth   =  3000;
    else
        mrsiReconParams.MaxSweepwidth   = 1500 ;
    end
else
    mrsiReconParams.MaxSweepwidth   = 1E12 ;
end


% LIPID SUPPRESSION PARAMS
mrsiReconParams.LipidRemovalParams.LipidRemStrength = LipidRemStrength;  % strength of the Lipid suppresion : 0 - 2 (unitless)
mrsiReconParams.LipidRemovalParams.LipidRemExpo = 1+ (LipidRemStrength>1)*(LipidRemStrength-1) ;
mrsiReconParams.LipidRemovalParams.TolLipRMDamping =0.92^LipidRemStrength; % tolerance for attenuation of overall spectral signal when applying the lipid suppresion (paramter will be determined this way)

mrsiReconParams.LipidMinPPM=-3.0; %ppm range used to identify lipid signal (in minus ppm )
mrsiReconParams.LipidMaxPPM=0.0;

% WATER SIGNal PROCESSING AND COIL SENSITIVITY PARAMS
mrsiReconParams.NbPtForWaterPhAmp= 3; %First pt of the time series taken for Amplitude and Phase calculation
mrsiReconParams.ESPIRIT_kernel=[5,5,5];
mrsiReconParams.SetInitFreqMapToZero = 0; %Set the initial B0 field map to 0 (discard the B0 fieldmap from the water data)
mrsiReconParams.LoadNiftiFreqMap = 0; %load B0 map from './B0MapMagn.nii.gz' and './B0MapPhase.nii.gz' (should be registered on MRSI resolution with FLIRT first);
mrsiReconParams.WaterRefFreqSearchRange =25*mrsiReconParams.FieldStrength;

mrsiReconParams.AutoMask_Quantile=0.9; % quantile used to compute the lipid and image mask. The higher the quantile the tighter the masks.
mrsiReconParams.GaussianSigma= 2.5; %Size in voxel of the smoothing 3D kernel for different fast reconstruction steps


% WATER RESIDUAL SUPPRESSION PARAMS

if mrsiReconParams.FieldStrength==7
    mrsiReconParams.FiltParam.Water_minFreq=-150;%hz
    mrsiReconParams.FiltParam.Water_maxFreq=150;%hz
else
    mrsiReconParams.FiltParam.Water_minFreq=-75;%hz
    mrsiReconParams.FiltParam.Water_maxFreq=75;%hz
end
mrsiReconParams.FiltParam.Comp=WaterSuppComp; % 0 for no Water residual suppression

mrsiReconParams.SideBandRMRank=0; % Rank for removal of the Water side band (assymetric residual Water signal peaks ) in the spectrum. Preleminary method. 0 for no Sideband removal. 64 is max
%mrsiReconParams.FiltParam.FilterFreq=0;%5;%hz


% Low-Rank TGV Recon PARAMS
mrsiReconParams.HammingConsistencyFilter=1; %Increase the recon quality for low SNR data
mrsiReconParams.LRTGVModelParams.maxit = 1500; % Number of iteration to perform for the LR-TGV recon
mrsiReconParams.LRTGVModelParams.check_it=100; %iteration step where iterative metrics are computed
mrsiReconParams.LRTGVModelParams.Plot_it=300;%;%iteration step where figures are made
mrsiReconParams.LRTGVModelParams.CorrB0Map_it=25; %iteration step where B0 map is corrected based on metabolite (MRSI) signal during recon (-999 disable this correction)
mrsiReconParams.LRTGVModelParams.CorrB0Map_Maxcount=30;% maximum number of  B0 map correction
mrsiReconParams.LRTGVModelParams.CorrB0Map_MaxPPM=-1.8;% ppm range for the B0 map correction
mrsiReconParams.LRTGVModelParams.CorrB0Map_MinPPM=-3.5;
mrsiReconParams.LRTGVModelParams.Orthogonalize_it=25;%iteration step where component are orthogonalized
mrsiReconParams.LRTGVModelParams.reduction=100; %Soften and accelerate connvergence by dynamically change the regularization during iterations
mrsiReconParams.LRTGVModelParams.Momentum=1./3; %Momentum given to the gradient during the descent. The higher the faster but the most unstable. Do not exceed 1/2.
mrsiReconParams.LRTGVModelParams.min_SpectStep=1/32;
mrsiReconParams.LRTGVModelParams.max_SpectStep=1/2;
mrsiReconParams.LRTGVModelParams.min_taup=1/32;
mrsiReconParams.LRTGVModelParams.max_taup=1/8;
mrsiReconParams.LRTGVModelParams.DualPrimalTauFact=2;

%Make directory for log file
if ~exist(mrsiReconParams.Log_Dir);mkdir(mrsiReconParams.Log_Dir);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and filter MRSI Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
fprintf('Reading MRSI  Data file. \n');
load(mrsiReconParams.MRSI_Filename);

mrsiReconParams.mrsiData_ctkkk = mrsiData_ctkkk;
mrsiReconParams.kmask = kmask;
mrsiReconParams.mrProt = mrProt;
clear mrsiData_ctkkk kmask mrProt

mrsiReconParams.mrProt.samplerate=1.0/mrsiReconParams.mrProt.DwellTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading TWIX file of the low resolution Water measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Reading Water Data. \n');
load(mrsiReconParams.Water_Filename);

mrsiReconParams.Water_ctkkk = Water_ctkkk;
mrsiReconParams.Water_kmask = Water_kmask;
mrsiReconParams.Water_mrProt = Water_mrProt;
clear Water_ctkkk Water_kmask Water_mrProt

mrsiReconParams.Water_mrProt.samplerate=1.0/mrsiReconParams.Water_mrProt.DwellTime;

if(numel(mrsiReconParams.mrsiData_ctkkk)<numel(mrsiReconParams.Water_ctkkk))
    warning('Water dataset contains more elements than the Main MRSI dataset! Files were probably swapped.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compress coils
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( mrsiReconParams.NbCoilCompression<999 ) & ( mrsiReconParams.NbCoilCompression < size(mrsiReconParams.Water_ctkkk,1) )

    K=mrsiReconParams.NbCoilCompression;
    [NbCoil Nt M N Slc ] = size(mrsiReconParams.Water_ctkkk);

    fprintf(makeSectionDisplayHeader(['Coil compression: going from ',num2str(NbCoil),' to ',num2str(K),' coils.']));

    Data_kc=permute(reshape(sqz(sum(mrsiReconParams.Water_ctkkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:,:),2)),[NbCoil M*N*Slc]),[2 1]);

    [~,S,V_cK]=svd(Data_kc,0);
    V_cK=single(V_cK(:,1:K));

    %Water data
    Temp=zeros([K Nt M N Slc ],class(mrsiReconParams.Water_ctkkk));
    for NewCoil=1:K
        Temp(NewCoil,:,:,:,:)=sum(single(mrsiReconParams.Water_ctkkk).*reshape(V_cK(:,NewCoil),[NbCoil 1 1 1 1]),1);
    end
    mrsiReconParams.Water_ctkkk = Temp;

    %Main Data

    [NbCoil Nt M N Slc ] = size(mrsiReconParams.mrsiData_ctkkk);
    Temp=zeros([K Nt M N Slc ],class(mrsiReconParams.mrsiData_ctkkk));
    for NewCoil=1:K
        Temp(NewCoil,:,:,:,:)=sum(single(mrsiReconParams.mrsiData_ctkkk).*reshape(V_cK(:,NewCoil),[NbCoil 1 1 1 1]),1);
    end
    mrsiReconParams.mrsiData_ctkkk = Temp;
    clear Temp

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduce FID Bandwidth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mrsiReconParams.mrProt.samplerate > mrsiReconParams.MaxSweepwidth
    mrsiReconParams.mrProt.VSize=2*round(0.5*mrsiReconParams.MaxSweepwidth/mrsiReconParams.mrProt.samplerate*mrsiReconParams.mrProt.VSize);
    mrsiReconParams.mrProt.samplerate = mrsiReconParams.MaxSweepwidth;
    mrsiReconParams.mrProt.DwellTime=1.0/mrsiReconParams.mrProt.samplerate;

    [NbCoil NbTime M N Slc] = size(mrsiReconParams.mrsiData_ctkkk);
    fprintf(['Reducing the number of points in the FID from ',num2str(NbTime),' to ',num2str(mrsiReconParams.mrProt.VSize),' ...\n']);
    Xd=1:(NbTime);
    Xm=linspace(1,NbTime,mrsiReconParams.mrProt.VSize);
    try
        mrsiReconParams.mrsiData_ctkkk=permute(interp1(Xd,permute(mrsiReconParams.mrsiData_ctkkk,[2,1,3,4,5]),Xm,'spline'),[2,1,3,4,5]);
    catch e

        mrsiData_ctkkk = zeros([NbCoil mrsiReconParams.mrProt.VSize M N Slc],class(mrsiReconParams.mrsiData_ctkkk));
        for c=1:size(mrsiReconParams.mrsiData_ctkkk,1)
            mrsiData_ctkkk(c,:,:,:,:)=interp1(Xd,sqz(mrsiReconParams.mrsiData_ctkkk(c,:,:,:,:)),Xm,'spline');
        end
        mrsiReconParams.mrsiData_ctkkk = mrsiData_ctkkk;
        clear mrsiData_tckkk
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute usefull variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDims = ndims(mrsiReconParams.mrsiData_ctkkk)-2;
Size_data=size(mrsiReconParams.mrsiData_ctkkk );

mrsiReconParams.spa_index = repmat({':'}, 1, nDims);
mrsiReconParams.ppm=(-4.7+((0:(mrsiReconParams.mrProt.VSize-1))*mrsiReconParams.mrProt.samplerate/(mrsiReconParams.mrProt.VSize*mrsiReconParams.mrProt.NMRFreq)));
[~,mrsiReconParams.MaxPPM_pt]=min(abs(mrsiReconParams.MaxPPM   - mrsiReconParams.ppm));
[~,mrsiReconParams.MinPPM_pt]=min(abs(mrsiReconParams.MinPPM   - mrsiReconParams.ppm));

mrsiReconParams.SlicesOversampling = round(mrsiReconParams.mrProt.Nslc*(mrsiReconParams.mrProt.FoV3D-mrsiReconParams.mrProt.VOI3D)/mrsiReconParams.mrProt.FoV3D);

fprintf( ['Slice oversampling estimated to be ',num2str( mrsiReconParams.SlicesOversampling),' slices.\n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do Water Supression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VisualizeKmasks( mrsiReconParams.mrsiData_ctkkk,mrsiReconParams.kmask, mrsiReconParams.Water_ctkkk ,mrsiReconParams.Water_kmask,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData]);
if mrsiReconParams.FiltParam.Comp>0
    fprintf(makeSectionDisplayHeader(['Starting Water Removal by HSVD with ',num2str(mrsiReconParams.FiltParam.Comp),' components ...']));
    Filtered_ctkkk=zeros(size(mrsiReconParams.mrsiData_ctkkk));
    parfor coil=1:Size_data(1);
        [Filtered_ctkkk(coil,:,:,:,:)] = WaterSuppression(squeeze(mrsiReconParams.mrsiData_ctkkk(coil,:,:,:,:)),mrsiReconParams.mrProt, mrsiReconParams.FiltParam, mrsiReconParams.kmask);
    end
    mrsiReconParams.mrsiData_ctkkk= Filtered_ctkkk;
    clear Filtered_ctkkk
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute other usefull variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Size_data=size(mrsiReconParams.mrsiData_ctkkk);

if mrsiReconParams.SpinEcho
    mrsiReconParams.NbMissingPoints=0;
    fprintf( ['SpinEcho sequence, Echo should be at first time point. No correction Needed.\n']);
else
    mrsiReconParams.NbMissingPoints=ceil(mrsiReconParams.mrProt.samplerate*mrsiReconParams.AcqDelay);%;
    fprintf( [num2str(mrsiReconParams.NbMissingPoints) ' calculated missing points.\n']);
end


if ~isfield(mrsiReconParams.mrProt,'Nlines')
    mrsiReconParams.mrProt.Ncol=size(mrsiReconParams.mrsiData_ctkkk,3);
    mrsiReconParams.mrProt.Nlines=size(mrsiReconParams.mrsiData_ctkkk,4);
    mrsiReconParams.mrProt.Nslc=size(mrsiReconParams.mrsiData_ctkkk,5);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the UnderSampling pattern in k-space if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mrsiReconParams.Originalkmask=mrsiReconParams.kmask;
if mrsiReconParams.UndersamplingF < 1
    [US_MASK]=Make_undersampled_mask(1.0,mrsiReconParams.UndersamplingF,mrsiReconParams.kmask,0.2);
    mrsiReconParams.kmask=mrsiReconParams.kmask.*US_MASK;

    s=[mrsiReconParams.Log_Dir,filesep,'Kmasks_Original_US_CS_',Data_name,'.ps'];
    if exist(s);delete(s);end
    figs=figure('visible','off');

    imagesc(Vol2Image(fftshift(mrsiReconParams.Originalkmask)));
    title('Original K-mask');
    print(figs, '-append', '-dpsc2', s);
    figs=figure('visible', 'off');
    imagesc(Vol2Image(fftshift(mrsiReconParams.kmask)));
    title('Random undersampled K-mask');
    print(figs, '-append', '-dpsc2', s);
    close all;

    fprintf(makeSectionDisplayHeader('Apply the Under-sampling pattern k-space mask ...'));
    kmask_ctkkk=permute(repmat(mrsiReconParams.kmask,[1,1,1,Size_data(1),Size_data(2)]),[4,5,1,2,3]);
    mrsiReconParams.mrsiData_ctkkk=mrsiReconParams.mrsiData_ctkkk.*kmask_ctkkk;

    clear kmask_ctkkk US_MASK
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Coil Sens Profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( 'Compute Anatomical masks ...\n');
[ mrsiReconParams.ImMask, mrsiReconParams.BrainMask,mrsiReconParams.SkMask ] = MakeMasks(mrsiReconParams );

fprintf( 'Compute Sensitivity profiles ...\n');
[mrsiReconParams.SENSE,mrsiReconParams.WaterFreqMap, mrsiReconParams.Water_trrr, mrsiReconParams.Water_rrr]=ComputeSENSEProfilesWithESPIRiT(mrsiReconParams);

VisualizeKmasks( mrsiReconParams.mrsiData_ctkkk,mrsiReconParams.kmask, mrsiReconParams.Water_ctkkk ,mrsiReconParams.Water_kmask,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData]);
VisualizeData( mrsiReconParams.mrsiData_ctkkk, mrsiReconParams.Water_ctkkk ,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData '_OriginalData']);
VisualizeMasks( mrsiReconParams.Water_ctkkk,mrsiReconParams.ImMask ,mrsiReconParams.BrainMask, mrsiReconParams.SkMask,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData '_WaterData']);
VisualizeCombinedData(mrsiReconParams.mrsiData_ctkkk, mrsiReconParams.Water_rrr ,mrsiReconParams,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData '_OriginalData']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute  Hann Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mrsiReconParams.HammingConsistencyFilter>0
    Ham=mrsiReconParams.HammingConsistencyFilter; %Haming filter strenght. 1= Max (Hann), 0=min (no filtering)
    fprintf( ['Compute and using Data-fidelity Hamming Filter with param',num2str(Ham) , ' to data and create Recon HannOperator...\n']);
    [X,Y,Z] = ndgrid(1:Size_data(3), 1:Size_data(4),1:Size_data(5));
    xc=floor(Size_data(3)/2)+1;yc=floor(Size_data(4)/2)+1;zc=floor(Size_data(5)/2)+1;
    temp = ((2*(X-xc)/Size_data(3)).^2 + (2*(Y-yc)/Size_data(4)).^2+ (2*(Z-zc)/Size_data(5)).^2).^0.5 ;

    mrsiReconParams.HKernel= fftshift( (1-Ham/2)+Ham/2*cos(pi*temp) );
    clear X Y Z xc xy xz temp
else
    fprintf( 'No  Data-fidelity Hamming Filter used ...\n');
    mrsiReconParams.HKernel=ones(Size_data(3),Size_data(4),Size_data(5));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Raw Data figures and check the Skull Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VisualizeData( mrsiReconParams.mrsiData_ctkkk, mrsiReconParams.Water_ctkkk ,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData '_BeforeLipidSupp']);

if (mrsiReconParams.LipidRemovalParams.LipidRemStrength>0)

    SENSE_c1rrr=reshape(mrsiReconParams.SENSE,[size(mrsiReconParams.SENSE,1) 1 size(mrsiReconParams.SENSE,2) size(mrsiReconParams.SENSE,3) size(mrsiReconParams.SENSE,4)]).*reshape(mrsiReconParams.ImMask,[1 1 size(mrsiReconParams.ImMask)]);
    Lipids_tkkk=squeeze(fft(fft(fft(sum(conj(SENSE_c1rrr).*ifft(ifft(ifft(mrsiReconParams.mrsiData_ctkkk,[],3),[],4),[],5),1),[],3),[],4),[],5));
    clear SENSE_c1rrr

    Lipid_Vol=squeeze(sum(abs(ifft(ifft(ifft(Lipids_tkkk,[],4),[],3),[],2).^2),1));

    Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));
    SkMask_temp=(Lipid_Vol>mean(Lipid_Vol(:)));%(Lipid_Vol>0.5*mean(Lipid_Vol(:)));
    overlap=(SkMask_temp.*mrsiReconParams.SkMask);
    diff=(SkMask_temp.*~mrsiReconParams.SkMask);

    mrsiReconParams.ResidWaterMask=diff.*mrsiReconParams.BrainMask;
    mrsiReconParams.SkMask=mrsiReconParams.SkMask+mrsiReconParams.ResidWaterMask;

    VisualizeMasks( mrsiReconParams.mrsiData_ctkkk, mrsiReconParams.ImMask ,mrsiReconParams.BrainMask, mrsiReconParams.SkMask,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData '_MRSIData_WaterResdiual-in-SkullMask' ]);
    if (sum(diff(:))/sum(overlap(:)))>2
        fprintf([ 'Difference / Overlap between Skull mask and Signal = ',num2str(sum(diff(:))/sum(overlap(:))), ' ( should be < 0.1 )\n']);
        error('There is too much signal in the Brain mask in comparison to the skull mask. There is obviously a mismatch between the masks and the MRSI data. Recon will not work. Stopping here!')
    elseif (sum(diff(:))/sum(overlap(:)))>0.05
        warning('A Lot of signal outside the skull mask before lipid suppression. There is probably a lot of unsuppressed water.')
        fprintf([ 'Difference / Overlap between Skull mask and Signal = ',num2str(sum(diff(:))/sum(overlap(:))), ' ( should be < 0.1 )\n']);
    end
end

mrsiReconParams=rmfield(mrsiReconParams, {'Water_ctkkk'}); % Not used during for the rest of the code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restrict MRSI data to the desired frequency window for the recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['Keeping points ',num2str(mrsiReconParams.MinPPM_pt),' to ',num2str(mrsiReconParams.MaxPPM_pt),' for reconstruction.\n']);

OriginalFullData_cfkkk=fft(mrsiReconParams.mrsiData_ctkkk ,[],2);
clear Time_1t111

mrsiReconParams=rmfield(mrsiReconParams, {'mrsiData_ctkkk'}); % Not used during for the rest of the code

if  (mrsiReconParams.LipidRemovalParams.LipidRemStrength>0)
    mrsiData_ckkkt=permute(ifft( OriginalFullData_cfkkk ,[],2), [1,3,4,5,2]);
else
    mrsiData_ckkkt=OriginalFullData_cfkkk(:,1:mrsiReconParams.MaxPPM_pt,:,:,:);
    mrsiData_ckkkt(:,1:mrsiReconParams.MinPPM_pt,:,:,:)=0;
    mrsiData_ckkkt=permute(ifft( mrsiData_ckkkt ,[],2), [1,3,4,5,2]);

    %Remove negative time signal
    fprintf("Filter out non-causal time signal before Recon...\n");
    Fs = mrsiReconParams.mrProt.samplerate*size(mrsiData_ckkkt,5)/mrsiReconParams.mrProt.VSize;

     mrsiData_ckkkt = permute(HSVD_TimeCausalityFilter(permute(mrsiData_ckkkt,[5,1,2,3,4]),Fs,min([64,floor(size(mrsiData_ckkkt,5)/4)])),[2,3,4,5,1]);

end

OriginalFullData_trrr=zeros(size(OriginalFullData_cfkkk,2),size(OriginalFullData_cfkkk,3),size(OriginalFullData_cfkkk,4),size(OriginalFullData_cfkkk,5));

[ ~,  M, N, O] = size(OriginalFullData_trrr);

SENSE_c1rrr=reshape(mrsiReconParams.SENSE,[ size(mrsiReconParams.SENSE,1) 1 M N O]);
OriginalFullData_trrr= conj(SENSE_c1rrr).*ifft(ifft(ifft(ifft(OriginalFullData_cfkkk,[],2),[],3),[],4),[],5);
OriginalFullData_trrr =squeeze(sum(OriginalFullData_trrr,1));

clear OriginalFullData_cfkkk SENSE_ctrrr;

if ~exist(mrsiReconParams.Results_Dir);mkdir(mrsiReconParams.Results_Dir);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOW RANK TGV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(makeSectionDisplayHeader('Low-Rank TGV reconstruction...'));

if mrsiReconParams.LipidRemovalParams.LipidRemStrength>0
    [mrsiReconParams.ReconDataShort_trrr, ~, ~ ,mrsiReconParams.U_rrrc, mrsiReconParams.V_tc, mrsiReconParams.S,mrsiReconParams.LipidRM_ff,mrsiReconParams.DynFreqMap] = LowRankTGV_ModelReconWithLipids(mrsiData_ckkkt,mrsiReconParams);
else
    [mrsiReconParams.ReconDataShort_trrr,mrsiReconParams.U_rrrc,mrsiReconParams.V_tc,mrsiReconParams.S,mrsiReconParams.DynFreqMap,mrsiReconParams.CostFunVal] = LowRankTGV_ModelRecon(mrsiData_ckkkt,mrsiReconParams);
end

clear TempRes_rrrt Recon_US  mrsiData_ckkkt;

if mrsiReconParams.LipidRemovalParams.LipidRemStrength>0
    OriginalSize_ReconData_frrr=fft(OriginalFullData_trrr,[],1); %frrr
    OriginalSize_ReconData_frrr(mrsiReconParams.MinPPM_pt:mrsiReconParams.MaxPPM_pt,:,:,:) = fft(mrsiReconParams.ReconDataShort_trrr,[],1);
else
    OriginalSize_ReconData_frrr=fft(OriginalFullData_trrr,[],1); %frrr
    OriginalSize_ReconData_frrr(1:mrsiReconParams.MaxPPM_pt,:,:,:)=fft(mrsiReconParams.ReconDataShort_trrr,[],1); %frrr
end


%Remove signal in time <0 
fprintf("Filter out non-causal time signal...\n");
Fs = mrsiReconParams.mrProt.samplerate*size(OriginalSize_ReconData_frrr,1)/mrsiReconParams.mrProt.VSize;
OriginalSize_ReconData_frrr = fft(HSVD_TimeCausalityFilter(ifft(OriginalSize_ReconData_frrr,[],1),Fs,min([64,floor(size(OriginalSize_ReconData_frrr,1)/4)])),[],1);

mrsiReconParams.Recon_OriginalSize_Data_trrr=ifft(OriginalSize_ReconData_frrr,[],1);%trrr
clear Time_trrr

clear OriginalSize_ReconData_frrr;

ReconData_frrr = fft(mrsiReconParams.Recon_OriginalSize_Data_trrr,[],1);

clear OriginalFullData_trrr

%%%%
%Add first missing points

Size_data=size(ReconData_frrr);
NaNData=[repmat(NaN,[mrsiReconParams.NbMissingPoints,Size_data(2)*Size_data(3)*Size_data(4)]); reshape(ifft(ReconData_frrr,[],1),Size_data(1),[]) ];
ReconData_frrr = fft(reshape(fillgaps(double(NaNData)),[ size(NaNData,1) Size_data(2:4)] ),[],1);

mrsiReconParams.Recon_Data_trrr=ifft( ReconData_frrr,[],1);
clear  ReconData_frrr

mrsiReconParams.mrProt.VSize=Size_data(1);
mrsiReconParams.ppm=(-4.7+((1:mrsiReconParams.mrProt.VSize)*mrsiReconParams.mrProt.samplerate/(mrsiReconParams.mrProt.VSize*mrsiReconParams.mrProt.NMRFreq)));

[~,mrsiReconParams.MaxPPM_pt]=min(abs(mrsiReconParams.MaxPPM - mrsiReconParams.ppm));
[~,mrsiReconParams.MinPPM_pt]=min(abs(mrsiReconParams.MinPPM - mrsiReconParams.ppm));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save  data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reconResults= mrsiReconParams;
clear mrsiReconParams

reconResults.MainRDA_Filename=[reconResults.Results_Dir,'/MRSI_CSSENSELR_Recon_',reconResults.NameData,'.rda'];
reconResults.OrigRDA_Filename=[reconResults.Results_Dir,'/MRSI_LipidandWaterRemoved_NoRecon_',reconResults.NameData,'.rda'];
reconResults.WaterRDA_Filename=[reconResults.Results_Dir,'/MRSI_Water_',reconResults.NameData,'.rda'];

reconResults.Recon_Data_trrr = single(reconResults.Recon_Data_trrr);
reconResults.Water_trrr = single(reconResults.Water_trrr);
reconResults.SENSE = single(reconResults.SENSE);
reconResults.ReconDataShort_trrr = single(reconResults.ReconDataShort_trrr);

if ~exist(reconResults.Results_Dir);mkdir(reconResults.Results_Dir);end

fprintf('Saving Data. \n');
save( [reconResults.Results_Dir,'/MRSI_MatlabReconResult_',reconResults.NameData,'.mat'], 'reconResults', '-v7.3');

end
