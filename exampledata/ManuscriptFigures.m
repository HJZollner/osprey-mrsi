% This script automatically the raw figures of the analysis of the Osprey MRSI example
% data which is presented in the manuscript 'Fully automated open-source
% analysis and interactive visualization of magnetic resonance spectroscopic imaging (MRSI) data in Osprey-MRSI'
% by Zollner et al. 2026.
% 
%   PREREQUISITS:
%   Osprey MRSI Matlab path (https://github.com/HJZollner/osprey-mrsi)
%   You have to run the RunAllDatasets.m script located in the 'osprey-mrsi/exampledata' folder first.
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2025-03-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2026-03-25: First version of the code.

assert(exist('OspreyMRSI')>0,'Please ensure that the osprey-mrsi folder (and subfolders) have been added to the MATLAB path')
assert(strlength(which((fullfile('exampledata','mrsi','Philips','TE_15_braino_phantom','derivatives','jobMRSI_TE_30_in_vitro_P_NII_ax.mat'))))>0, 'Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
%% In Vitro Scan Philips
assert(strlength(which((fullfile('exampledata','mrsi','Philips','TE_15_braino_phantom','derivatives','jobMRSI_TE_30_in_vitro_P_NII_ax.mat'))))>0, 'Philips in vitro results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','Philips','TE_15_braino_phantom','derivatives','jobMRSI_TE_15_in_vitro_SPARSDAT.mat')));

out = osp_plotQuickmapsOverlay(MRSCont,'w','H2O',1,5,1,0,'T1w_rMRSI',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_Overlay_QuickWater.pdf','pdf');
saveas(gcf,'SI_7_Overlay_QuickWater.fig','fig');
close(gcf);

out = osp_plotSpecAndLocMRSI(MRSCont,[11,15,1;], 'T1w_rMRSI','OspreyProcess','w','Fit3D',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_Voxel_1.pdf','pdf');
saveas(gcf,'SI_7_Voxel_1.fig','fig');
close(gcf);

out = osp_plotSpecAndLocMRSI(MRSCont,[11,15,2;], 'T1w_rMRSI','OspreyProcess','w','Fit3D',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_Voxel_2.pdf','pdf');
saveas(gcf,'SI_7_Voxel_2.fig','fig');
close(gcf);

out = osp_plotSpecAndLocMRSI(MRSCont,[11,15,3;], 'T1w_rMRSI','OspreyProcess','w','Fit3D',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_Voxel_3.pdf','pdf');
saveas(gcf,'SI_7_Voxel_3.fig','fig');
close(gcf);

out = osp_plotSpecAndLocMRSI(MRSCont,[11,15,4;], 'T1w_rMRSI','OspreyProcess','w','Fit3D',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_Voxel_4.pdf','pdf');
saveas(gcf,'SI_7_Voxel_4.fig','fig');
close(gcf);

out = osp_plotSpecAndLocMRSI(MRSCont,[11,15,5], 'T1w_rMRSI','OspreyProcess','w','Fit3D',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_Voxel_5.pdf','pdf');
saveas(gcf,'SI_7_Voxel_5.fig','fig');
close(gcf);

out = osp_plotSpecAndLocMRSI(MRSCont,[11,15,1;11,15,2;11,15,3;11,15,4;11,15,5], 'T1w_rMRSI','OspreyProcess','w','Fit3D',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_Voxel_all.pdf','pdf');
saveas(gcf,'SI_7_Voxel_all.fig','fig');
close(gcf);

%% Phantom scan GE ACR
assert(strlength(which((fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','derivatives_sag','jobMRSI_TE_30_in_vitro_P_NII_sag.mat'))))>0, 'GE saggital in vitro results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','derivatives_sag','jobMRSI_TE_30_in_vitro_P_NII_sag.mat')));
out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI', 1, 0, 0, 0,1,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_GE_ACR_mask_sag.pdf','pdf');
saveas(gcf,'SI_7_GE_ACR_mask_sag.fig','fig');
close(gcf);

assert(strlength(which((fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','derivatives_ax','jobMRSI_TE_30_in_vitro_P_NII_ax.mat'))))>0, 'GE axial in vitro results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','derivatives_ax','jobMRSI_TE_30_in_vitro_P_NII_ax.mat')));
out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI', 1, 0, 0, 0,1,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_GE_ACR_mask_ax.pdf','pdf');
saveas(gcf,'SI_7_GE_ACR_mask_ax.fig','fig');
close(gcf);

%% In Vitro scan Siemens ACR
assert(strlength(which((fullfile('exampledata','mrsi','Siemens','TE_40_phantom_acr','derivatives','jobMRSI_TE_40_in_vitro_RDA_NII.mat'))))>0, 'Siemens in vitro results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','Siemens','TE_40_phantom_acr','derivatives','jobMRSI_TE_40_in_vitro_RDA_NII.mat')))

out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI', 1, 0, 0, 0,2,2);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_7_Siemens_ACR_mask_ax.pdf','pdf');
saveas(gcf,'SI_7_Siemens_ACR_mask_ax.fig','fig');
close(gcf);

%% In vivo medium-TE MRSI dataset figures
assert(strlength(which((fullfile('exampledata','mrsi','Philips','TE_70','derivatives','jobMRSI_TE_70_SPARSDAT.mat'))))>0, 'Philips in vitro medium-TE results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','Philips','TE_70','derivatives','jobMRSI_TE_70_SPARSDAT.mat')));

% Segmentation & Coregistration
slice = 4;
out = osp_plotSegmentationOverlay(MRSCont,'fCSF',slice,slice,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_fCSF.pdf','pdf');
saveas(gcf,'Figure_2_fCSF.fig','fig');
close(gcf);
out = osp_plotSegmentationOverlay(MRSCont,'fWM',slice,slice,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_fWM.pdf','pdf');
saveas(gcf,'Figure_2_fWM.fig','fig');
close(gcf);
out = osp_plotSegmentationOverlay(MRSCont,'fGM',slice,slice,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_fGM.pdf','pdf');
saveas(gcf,'Figure_2_fGM.fig','fig');
close(gcf);
out = osp_plotSegmentationOverlay(MRSCont,'fGM',slice,slice,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_seg_colorbar.pdf','pdf');
saveas(gcf,'Figure_2_seg_colorbar.fig','fig');
close(gcf);
out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSIloc', 0, 0, 1, 1,slice,slice);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_mask.pdf','pdf');
saveas(gcf,'Figure_2_mask.fig','fig');
close(gcf);

% Spectra [25,23,4;15,19,4]
out = osp_plotSpecAndLocMRSI(MRSCont,[25,23,4;15,19,4], 'T1w_rMRSIloc','OspreyLoad','AFID','Fit3D',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_Raw.pdf','pdf');
saveas(gcf,'Figure_2_Raw.fig','fig');
close(gcf);

out = osp_plotSpecAndLocMRSI(MRSCont,[25,23,4;15,19,4], 'T1w_rMRSIloc','OspreyProcess','AFID','Fit3D',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_Proc.pdf','pdf');
saveas(gcf,'Figure_2_Proc.fig','fig');
close(gcf);

out = osp_plotSpecAndLocMRSI(MRSCont,[25,23,4;15,19,4], 'T1w_rMRSIloc','OspreyFit','metab','Fit1DStack',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_Fit.pdf','pdf');
saveas(gcf,'Figure_2_Fit.fig','fig');
close(gcf);

% QC maps
out = osp_plotMetabolitemapsOverlay(MRSCont,'GlobalQC','QC',slice,slice,1,1,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_GlobalQC_colorbar.pdf','pdf');
saveas(gcf,'Figure_2_GlobalQC_colorbar.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'GlobalQC','QC',slice,slice,1,1,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_GlobalQC.pdf','pdf');
saveas(gcf,'Figure_2_GlobalQC.fig','fig');
close(gcf);

% Metabolite maps
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tNAA_Acetyl',slice,slice,1,30,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tNAA.pdf','pdf');
saveas(gcf,'Figure_2_tNAA.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tNAA_Acetyl',slice,slice,1,30,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tNAA_cbar.pdf','pdf');
saveas(gcf,'Figure_2_tNAA_cbar.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tCr_methyl',slice,slice,1,16,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCr.pdf','pdf');
saveas(gcf,'Figure_2_tCr.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tCr_methyl',slice,slice,1,16,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCr_cbar.pdf','pdf');
saveas(gcf,'Figure_2_tCr_cbar.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tCho_methyl',slice,slice,1,4,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCho.pdf','pdf');
saveas(gcf,'Figure_2_tCho.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tCho_methyl',slice,slice,1,4,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCho_cbar.pdf','pdf');
saveas(gcf,'Figure_2_tCho_cbar.fig','fig');
close(gcf);

% CRLB
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tNAA_Acetyl',slice,slice,1,15,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tNAA_CRLBs.pdf','pdf');
saveas(gcf,'Figure_2_tNAA_CRLBs.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tNAA_Acetyl',slice,slice,1,15,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tNAA_cbar_CRLBs.pdf','pdf');
saveas(gcf,'Figure_2_tNAA_cbar_CRLBs.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tCr_methyl',slice,slice,1,15,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCr_CRLBs.pdf','pdf');
saveas(gcf,'Figure_2_tCr_CRLBs.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tCr_methyl',slice,slice,1,15,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCr_cbar_CRLBs.pdf','pdf');
saveas(gcf,'Figure_2_tCr_cbar_CRLBs.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tCho_methyl',slice,slice,1,15,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCho_CRLBs.pdf','pdf');
saveas(gcf,'Figure_2_tCho_CRLBs.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tCho_methyl',slice,slice,1,50,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCho_cbar_CRLBs.pdf','pdf');
saveas(gcf,'Figure_2_tCho_cbar_CRLBs.fig','fig');
close(gcf);

% Atlas
% Run the follwoing coment and pick thalamus at 33%
% out = osp_plotInteractiveAtlasAnalysis(MRSCont,1,'TissCorrWaterScaled',MRSCont.opts.MRSI.atlas.metabolites,'T1w_rMRSI')
% set(gcf,'Renderer', 'painters');
% saveas(gcf,'Figure_2_Atlas.pdf','pdf');
% close(gcf);


% Global
out = osp_plotGlobalConcentration(MRSCont,'tCho_methyl','TissCorrWaterScaled');
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_2_tCho_GlobalConc.pdf','pdf');
saveas(gcf,'Figure_2_tCho_GlobalConc.fig','fig');
close(gcf);

%% In vivo short-TE MRSI dataset figures.
assert(strlength(which((fullfile('exampledata','mrsi','Philips','TE_15','derivatives','jobMRSI_TE_15_SPARSDAT.mat'))))>0, 'Philips in vivo short-TE results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','Philips','TE_15','derivatives','jobMRSI_TE_15_SPARSDAT.mat')));

% Voxels Coreg
out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSIloc', 1, 0, 0, 0,1,5);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_mask.pdf','pdf');
saveas(gcf,'Figure_3_mask.fig','fig');
close(gcf);


out = osp_plotSpecAndLocMRSI(MRSCont,[12,11,2;15,13,2;13,16,2;12,18,2], 'T1w_rMRSIloc','OspreyFit','metab','Fit1DStack',0,0,0,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_Fit.pdf','pdf');
saveas(gcf,'Figure_3_Fit.fig','fig');
close(gcf);


% Metabolite maps
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tNAA',3,4,1,35,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tNAA.pdf','pdf');
saveas(gcf,'Figure_3_tNAA.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tNAA',3,4,1,35,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tNAA_cbar.pdf','pdf');
saveas(gcf,'Figure_3_tNAA_cbar.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tCr',3,4,1,20,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tCr.pdf','pdf');
saveas(gcf,'Figure_3_tCr.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tCr',3,4,1,20,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tCr_cbar.pdf','pdf');
saveas(gcf,'Figure_3_tCr_cbar.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tCho',3,4,1,6,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tCho.pdf','pdf');
saveas(gcf,'Figure_3_tCho.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','tCho',3,4,1,6,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tCho_cbar.pdf','pdf');
saveas(gcf,'Figure_3_tCho_cbar.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','mI',3,4,1,15,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_mI.pdf','pdf');
saveas(gcf,'Figure_3_mI.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','mI',3,4,1,15,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_mI_cbar.pdf','pdf');
saveas(gcf,'Figure_3_mI_cbar.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','Glx',3,4,1,25,0,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_Glx.pdf','pdf');
saveas(gcf,'Figure_3_Glx.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'TissCorrWaterScaled','Glx',3,4,1,25,0,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_Glx_cbar.pdf','pdf');
saveas(gcf,'Figure_3_Glx_cbar.fig','fig');
close(gcf);

% CRLB maps
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tNAA',3,4,1,15,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tNAA_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_tNAA_CRLBs.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tNAA',3,4,1,15,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tNAA_cbar_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_tNAA_cbar_CRLBs.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tCr',3,4,1,15,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tCr_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_tCr_CRLBs.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tCr',3,4,1,15,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tCr_cbar_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_tCr_cbar_CRLBs.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tCho',3,4,1,40,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tCho_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_tCho_CRLBs.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','tCho',3,4,1,40,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_tCho_cbar_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_tCho_cbar_CRLBs.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','mI',3,4,1,40,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_mI_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_mI_CRLBs.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','mI',3,4,1,40,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_mI_cbar_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_mI_cbar_CRLBs.fig','fig');
close(gcf);

out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','Glx',3,4,1,40,1,'T1w_rMRSIloc',0.8,0);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_Glx_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_Glx_CRLBs.fig','fig');
close(gcf);
out = osp_plotMetabolitemapsOverlay(MRSCont,'CRLBs','Glx',3,4,1,40,1,'T1w_rMRSIloc',0.8,1);
set(gcf,'Renderer', 'painters');
saveas(gcf,'Figure_3_Glx_cbar_CRLBs.pdf','pdf');
saveas(gcf,'Figure_3_Glx_cbar_CRLBs.fig','fig');
close(gcf);

% Atlas
% Run the follwoing coment and pick thalamus at 25%
% out = osp_plotInteractiveAtlasAnalysis(MRSCont,1,'TissCorrWaterScaled',MRSCont.opts.MRSI.atlas.metabolites,'T1w_rMRSI')
% set(gcf,'Renderer', 'painters');
% saveas(gcf,'Figure_3_Atlas.pdf','pdf');
% close(gcf);

%% In vivo GE MRSI datasets Supplementary Material
% short-TE MRSI dataset 30 ms
assert(strlength(which((fullfile('exampledata','mrsi','GE','TE_30','derivatives','jobMRSI_TE_30_P_NII.mat'))))>0, 'GE in vivo short-TE results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','GE','TE_30','derivatives','jobMRSI_TE_30_P_NII.mat')));
out = osp_plotSpecAndLocMRSI(MRSCont,[8,7,1;9,7,1;8,8,1;], 'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,2,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_8_Fit_TE30.pdf','pdf');
saveas(gcf,'SI_8_Fit_TE30.fig','fig');
close(gcf);

% medium-TE MRSI dataset 70 ms
assert(strlength(which((fullfile('exampledata','mrsi','GE','TE_70','derivatives','jobMRSI_TE_70_P_NII.mat'))))>0, 'GE in vivo medium-TE results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','GE','TE_70','derivatives','jobMRSI_TE_70_P_NII.mat')));
out = osp_plotSpecAndLocMRSI(MRSCont,[8,7,1;9,7,1;8,8,1;], 'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,2,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_8_Fit_TE70.pdf','pdf');
saveas(gcf,'SI_8_Fit_TE70.fig','fig');
close(gcf);

%% In vivo Siemens MRSI datasets Supplementary Material

% short-TE MRSI dataset 30 ms
assert(strlength(which((fullfile('exampledata','mrsi','Siemens','TE_35','derivatives','jobMRSI_TE_35_RDA_NII.mat'))))>0, 'Siemens in vivo short-TE results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','Siemens','TE_35','derivatives','jobMRSI_TE_35_RDA_NII.mat')));
out = osp_plotSpecAndLocMRSI(MRSCont,[9,8,1;10,8,1;11,8,1;], 'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,2,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_9_Fit_TE35_1.pdf','pdf');
saveas(gcf,'SI_9_Fit_TE35_1.fig','fig');
close(gcf);
out = osp_plotSpecAndLocMRSI(MRSCont,[9,9,1;10,9,1;11,9,1;], 'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,2,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_9_Fit_TE35_2.pdf','pdf');
saveas(gcf,'SI_9_Fit_TE35_2.fig','fig');
close(gcf);
out = osp_plotSpecAndLocMRSI(MRSCont,[9,8,1;10,8,1;11,8,1;9,9,1;10,9,1;11,9,1;], 'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,2,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_9_Fit_TE35_3.pdf','pdf');
saveas(gcf,'SI_9_Fit_TE35_3.fig','fig');
close(gcf);


% long-TE MRSI dataset TE 280 ms
assert(strlength(which((fullfile('exampledata','mrsi','Siemens','TE_280','derivatives','jobMRSI_TE_280_RDA_NII.mat'))))>0, 'Siemens in vivo long-TE results are missing. Please run the "osprey-mrsi/exampledata/RunAllDatasets.m" script before proceeding with the plots.')
load(which(fullfile('exampledata','mrsi','Siemens','TE_280','derivatives','jobMRSI_TE_280_RDA_NII.mat')));

out = osp_plotSpecAndLocMRSI(MRSCont,[9,8,1;10,8,1;11,8,1;], 'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,2,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_9_Fit_TE280_1.pdf','pdf');
saveas(gcf,'SI_9_Fit_TE280_1.fig','fig');
close(gcf);
out = osp_plotSpecAndLocMRSI(MRSCont,[9,9,1;10,9,1;11,9,1;], 'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,2,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_9_Fit_TE280_2.pdf','pdf');
saveas(gcf,'SI_9_Fit_TE280_2.fig','fig');
close(gcf);
out = osp_plotSpecAndLocMRSI(MRSCont,[9,8,1;10,8,1;11,8,1;9,9,1;10,9,1;11,9,1;], 'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,2,0,1);
set(gcf, 'units','normalized','outerposition',[0 0 1 0.66]);
set(gcf,'Renderer', 'painters');
saveas(gcf,'SI_9_Fit_TE280_3.pdf','pdf');
saveas(gcf,'SI_9_Fit_TE280_3.fig','fig');
close(gcf);
