% This script automatically runs the analysis of the Osprey MRSI example
% data which is presented in the manuscript 'Fully automated open-source
% analysis of magnetic resonance spectroscopic imaging (MRSI) data in Osprey'
% by Zollner et al. 2026.
% 
%   PREREQUISITS:
%   Osprey MRSI from the MRSI GitHub branch in Matlab path
%   (https://github.com/schorschinho/osprey/tree/MRSI)
%   All toolboxes required to run Osprey
%   (https://schorschinho.github.io/osprey/getting-started.html#system-requirements)
%   
%   Matlab Parallel Computing Toolbox to accelerate the modeling. If not available 
%   you have to change opts.MRSI.LCM.parallelComputing = 1; to 0 in all job
%   files.
%
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2025-03-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2025-03-06: First version of the code.


%% In Vitro datasets
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','Philips','TE_15_braino_phantom','jobMRSI_TE_15_in_vitro_SPARSDAT.m')),'11',1);
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','jobMRSI_TE_30_in_vitro_P_NII_ax.m')),'11',1);
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','jobMRSI_TE_30_in_vitro_P_NII_sag.m')),'11',1);
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','Siemens','TE_40_phantom_acr','jobMRSI_TE_40_in_vitro_RDA_NII.m')),'11',1);

%% In Vivo datasets
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','Philips','TE_70','jobMRSI_TE_70_SPARSDAT.m')),'11');
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','GE','TE_70','jobMRSI_TE_70_P_NII.m')),'11');
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','Siemens','TE_280','jobMRSI_TE_280_RDA_NII.m')),'11');
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','Philips','TE_15','jobMRSI_TE_15_SPARSDAT.m')),'11');
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','GE','TE_30','jobMRSI_TE_30_P_NII.m')),'11');
MRSCont = OspreyMRSI(which(fullfile('exampledata','mrsi','Siemens','TE_35','jobMRSI_TE_35_RDA_NII.m')),'11');



