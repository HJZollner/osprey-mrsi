function [MRSCont] = OspreyMRSIHTMLReport(MRSCont,kk,limits)
%% [MRSCont] = OspreyHTMLReport(MRSCont,kk)
%   This function creates a short HTML report of the processing and modeling
%   and should be called at the end of the analysis. It uses plotly to make
%   the results interactive.
%
%   USAGE:
%       MRSCont = OspreyHTMLReport(MRSCont,kk);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       kk          = subject index
%       limits      = ppm range of data
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-10-31)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-10-31: First version of the code.
%% Prepartion
if nargin < 3
    limits = [0,6];
end
varStruct.files_quickMaps_raw = [];
varStruct.names_quickMaps_raw = [];
varStruct.files_quickMaps_proc = [];
varStruct.names_quickMaps_proc = [];
varStruct.files_Quantification_Maps=[];
varStruct.names_Quantification_Maps=[];
varStruct.files_Quantification_Maps_QC=[];
varStruct.names_Quantification_Maps_QC=[];
varStruct.files_quickMaps_raw = [];
varStruct.names_quickMaps_raw = [];
varStruct.files_Quantification_Maps_QCfilt=[];
varStruct.names_Quantification_Maps_QCfilt=[];
varStruct.files_CRLB_Maps=[];
varStruct.names_CRLB_Maps=[];
varStruct.files_GlobalConc=[];
varStruct.names_GlobalConc=[];
varStruct.files_atlas=[];
varStruct.names_atlas=[];
varStruct.PosLoadSpec = [];
varStruct.CoregPos = [];
varStruct.SegPos = [];
varStruct.CoregSegPos = [];
varStruct.PosProcSpec = [];
varStruct.PosFitSpec = [];
varStruct.PosGlobalConc = [];
varStruct.PosAtlas = [];
varStruct.TargetSpec = [];

progressMsg('Starting OspreyMRSIHTMLReport for kk = %d', kk);
progressMsg('Limits: [%.3f %.3f]', limits(1), limits(2));

% Get colormaps and setup the inital path
progressMsg('Setting up output folders and subject information...');
colormaps = MRSCont.colormap;
ppmmin = 0.2;
ppmmax=4.2;

outputFolder    = fullfile(MRSCont.outputFolder,'Reports');
split_subject_path = strsplit(fileparts(MRSCont.files{kk}),filesep);
str_ind_sub = find(contains(split_subject_path,'sub'));
if ~isempty(str_ind_sub)
    sub_str = split_subject_path(str_ind_sub(1));
    sub_str = sub_str{1};
else
    sub_str = ['sub-' num2str(kk)];
end
outputFigures   = fullfile(MRSCont.outputFolder,'Reports','reportFigures',sub_str);
[foldername,filename,~]  = fileparts(MRSCont.files{kk});

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
if ~exist(outputFigures,'dir')
    mkdir(outputFigures);
end


[varStruct] = createInteractiveFigure(MRSCont,'CheckPlotly',varStruct,limits,outputFigures);


%% OspreyLoad
if MRSCont.flags.didLoadData

    progressMsg('Entering OspreyLoad section...');
    [varStruct] = createInteractiveFigure(MRSCont,'OspreyLoadSpectra',varStruct,limits,outputFigures);

    [varStruct] = createInteractiveFigure(MRSCont,'OspreyLoadQuickMaps',varStruct,limits,outputFigures);
      
    progressMsg('Finished OspreyLoad section.');
end

%% OspreyCoreg
if MRSCont.flags.didCoreg
    progressMsg('Entering OspreyCoreg section...');
       
    [varStruct] = createInteractiveFigure(MRSCont,'OspreyCoreg',varStruct,limits,outputFigures);

    progressMsg('Finished OspreyCoreg section.');
end

%% OspreySeg
if MRSCont.flags.didSeg
    progressMsg('Entering OspreySeg section...');
    
    [varStruct] = createInteractiveFigure(MRSCont,'OspreySeg',varStruct,limits,outputFigures);

    progressMsg('Finished OspreySeg section.');
end

%% Osprey Process
if MRSCont.flags.didProcess
    progressMsg('Entering OspreyProcess section...');

   [varStruct] = createInteractiveFigure(MRSCont,'OspreyProcessSpectra',varStruct,limits,outputFigures);
    
    progressMsg('OspreyProcess: preparing processed quick maps...');
    % Quick Maps Process
    field_names = fieldnames(MRSCont.opts.MRSI.quickMapsList{2});
    for ff = 1 : length(field_names)
        MRSCont.opts.MRSI.quickMaps.(field_names{ff}) = MRSCont.opts.MRSI.quickMapsList{2}.(field_names{ff});
    end

    
    [varStruct] = createInteractiveFigure(MRSCont,'OspreyProcessQuickMaps',varStruct,limits,outputFigures);

    [varStruct] = createInteractiveFigure(MRSCont,'OspreyProcessSNRMap',varStruct,limits,outputFigures);

    [varStruct] = createInteractiveFigure(MRSCont,'OspreyProcessFWHMMap',varStruct,limits,outputFigures);

    [varStruct] = createInteractiveFigure(MRSCont,'OspreyProcessGlobalQCMap',varStruct,limits,outputFigures);
    
end

%% OspreyFit
% Plot spectra
if MRSCont.flags.didFit
    progressMsg('Entering OspreyFit section...');
    
    [varStruct] = createInteractiveFigure(MRSCont,'OspreyFit',varStruct,limits,outputFigures);

    progressMsg('Finished OspreyFit section.');
end

%% OspreyQuantify
if MRSCont.flags.didQuantify
    progressMsg('Entering OspreyQuantify section...');

    [varStruct] = createInteractiveFigure(MRSCont,'OspreyQuantifyMap',varStruct,limits,outputFigures);

    [varStruct] = createInteractiveFigure(MRSCont,'OspreyQuantifyMapQC',varStruct,limits,outputFigures);

    [varStruct] = createInteractiveFigure(MRSCont,'OspreyQuantifyMapQCFilt',varStruct,limits,outputFigures);

    % CRLB maps
    progressMsg('OspreyQuantify: exporting CRLB maps...');
    [varStruct] = createInteractiveFigure(MRSCont,'OspreySeg',varStruct,limits,outputFigures);
    progressMsg('Finished OspreyQuantify section.');
end

%% OspreyOverview
if MRSCont.flags.didOverview
    progressMsg('Entering OspreyOverview section...');

    % Global Concentration Results
    if isfield(MRSCont.opts.MRSI,'OspreyQuantifyMapCRLB')
        progressMsg('OspreyOverview: exporting global concentration results...');
        [varStruct] = createInteractiveFigure(MRSCont,'OspreyOverviewGlobalConc',varStruct,limits,outputFigures);
    end

    % Atlas analysis results
    if isfield(MRSCont.opts.MRSI,'atlas')
        progressMsg('OspreyOverview: exporting atlas results...');
        [varStruct] = createInteractiveFigure(MRSCont,'OspreyOverviewAtlas',varStruct,limits,outputFigures);
    end
    progressMsg('Finished OspreyOverview section.');
end

%% Write report in HTML files
progressMsg('Starting final HTML report writing...');

logoPath=which(fullfile('graphics','osprey.png'));
copyfile(logoPath,fullfile(outputFigures, 'osprey.png'));

%Write as relative path
outputFigures   = fullfile('reportFigures',sub_str);

%write an html report:
progressMsg('Opening report file: %s', fullfile(outputFolder,[sub_str,'-report.html']));
fid=fopen(fullfile(outputFolder,[sub_str,'-report.html']),'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
fprintf(fid,'\n<body>');


if ~isempty(logoPath)
    fprintf(fid,'\n<img src= " %s " width="35" height="28"> <b> \tOsprey MRSI Analysis Report</b> ',fullfile(outputFigures, 'osprey.png'));
else
    fprintf(fid,'\n<b> Osprey MRSI Analysis Report</b>');
end
fprintf(fid,'\n<p><b>DATE:</b> %s \t <b>FILENAME:</b> %s </p>',date,filename);

if MRSCont.flags.didLoadData
    progressMsg('Writing OspreyLoad section to HTML report...');

    % OspreyLoad
    fprintf(fid,'\n<h2> Osprey Load</h2>');
    fprintf(fid,'\n<h3> Example Raw Spectra </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%.0fx" frameborder="0"></iframe>',fullfile(outputFigures, 'Raw.html'),varStruct.PosLoadSpec(4)*1.5);

    fprintf(fid,'\n<h3> Quickmaps Raw Data (Amplitude Integration)</h3>');
    for ff = 1 : length(varStruct.files_quickMaps_raw)
        fprintf(fid,'\n %s ', varStruct.names_quickMaps_raw{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',varStruct.files_quickMaps_raw{ff});
    end
end

if MRSCont.flags.didCoreg
    progressMsg('Writing OspreyCoreg section to HTML report...');

    % OspreyCoreg
    fprintf(fid,'\n<h2> Osprey Coregistration</h2>');
    fprintf(fid,'\n<h3> MRSI slice localization </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%.0fpx" frameborder="0"></iframe>',fullfile(outputFigures, 'Coreg.png'),varStruct.CoregPos(4)*1.8);
end

if MRSCont.flags.didSeg
    progressMsg('Writing OspreySeg section to HTML report...');

    % OspreySeg
    fprintf(fid,'\n<h2> Osprey Segmentation</h2>');
    fprintf(fid,'\n<h3> Outer mask + automated brain mask </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%.0fpx" frameborder="0"></iframe>',fullfile(outputFigures, 'CoregSeg.png'),varStruct.CoregSegPos(4)*1.8);

    fprintf(fid,'\n<h3> Tissue fraction maps + automated masks </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%.0fpx" frameborder="0"></iframe>',fullfile(outputFigures, 'Seg.html'),varStruct.PosAtlas(4)*1.5);
end

if MRSCont.flags.didProcess
    progressMsg('Writing OspreyProcess section to HTML report...');

    % OspreyProcess
    fprintf(fid,'\n<h2> Osprey Process</h2>');
    fprintf(fid,'\n<h3> Example Processed Spectra %s</h3>',varStruct.TargetSpec);
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%.0fx" frameborder="0"></iframe>',fullfile(outputFigures, 'Process.html'),varStruct.PosProcSpec(4)*1.5);

    fprintf(fid,'\n<h3> Quickmaps Processed Data (Amplitude Integration)</h3>');
    for ff = 1 : length(varStruct.files_quickMaps_proc)
        if ff ==  length(varStruct.files_quickMaps_proc)-2
            fprintf(fid,'\n<h3> Quality Metric Maps </h3>');
        end
        fprintf(fid,'\n %s ', varStruct.names_quickMaps_proc{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',varStruct.files_quickMaps_proc{ff});
    end
end

if MRSCont.flags.didFit
    progressMsg('Writing OspreyFit section to HTML report...');

    % OspreyFit
    fprintf(fid,'\n<h2> Osprey Fit</h2>');
    fprintf(fid,'\n<h3> Example Fits </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%.0fx" frameborder="0"></iframe>',fullfile(outputFigures, 'Fit.html'),varStruct.PosFitSpec(4)*1.5);
end

if MRSCont.flags.didQuantify
    progressMsg('Writing OspreyQuantify section to HTML report...');

    % OspreyQuantify
    fprintf(fid,'\n<h2> Osprey Quantify</h2>');
    fprintf(fid,'\n<h3> Metabolite maps (no QC applied) </h3>');
    for ff = 1 : length(varStruct.files_Quantification_Maps)
        fprintf(fid,'\n %s ', varStruct.names_Quantification_Maps{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',varStruct.files_Quantification_Maps{ff});
    end

    fprintf(fid,'\n<h3> Relative CRLB maps </h3>');
    for ff = 1 : length(varStruct.files_CRLB_Maps)
        fprintf(fid,'\n %s ', varStruct.names_CRLB_Maps{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',varStruct.files_CRLB_Maps{ff});
    end

    fprintf(fid,'\n<h3> QC filter maps </h3>');
    for ff = 1 : length(varStruct.files_Quantification_Maps_QC)
        fprintf(fid,'\n %s ', varStruct.names_Quantification_Maps_QC{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',varStruct.files_Quantification_Maps_QC{ff});
    end

    fprintf(fid,'\n<h3> Metabolite maps (QC filtered) </h3>');
    for ff = 1 : length(varStruct.files_Quantification_Maps_QCfilt)
        fprintf(fid,'\n %s ', varStruct.names_Quantification_Maps_QCfilt{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',varStruct.files_Quantification_Maps_QCfilt{ff});
    end
end

if MRSCont.flags.didOverview
    progressMsg('Writing OspreyOverview section to HTML report...');

    % OspreyOverview
    fprintf(fid,'\n<h2> Osprey Overview</h2>');
    fprintf(fid,'\n<h3> Atlas Analysis </h3>');
    for ff = 1 : length(varStruct.files_atlas)
        fprintf(fid,'\n <h4> %s </h4>', varStruct.names_atlas{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%.0fx" frameborder="0"></iframe>',varStruct.files_atlas{ff},varStruct.PosAtlas(4)*1.5);
    end

    fprintf(fid,'\n<h3> Global Concentrations </h3>');
    for ff = 1 : length(varStruct.files_GlobalConc)
        fprintf(fid,'\n <h4> %s </h4>', varStruct.names_GlobalConc{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%.0fx" frameborder="0"></iframe>',varStruct.files_GlobalConc{ff},varStruct.PosGlobalConc(4)*1.5);
    end
end

fprintf(fid,'\n</body>');
fprintf(fid,'\n</html>');
fclose(fid);

progressMsg('Finished writing HTML report.');
progressMsg('Completed OspreyMRSIHTMLReport for subject %s.', sub_str);
end


function [p] = cleanup_montages(p)
    progressMsg('cleanup_montages: starting...');

    p.layout.height = 200*2;
    p.layout.width = 840*2;
    p.layout.yaxis1.scaleanchor = 'x';
    p.layout.yaxis1.scaleratio = 1;
    p.layout.yaxis1.autorange = 'reversed';
    p.layout.xaxis1.showgrid = false;
    p.layout.yaxis1.showgrid = false;
    p.layout.xaxis1.zeroline = false;
    p.layout.yaxis1.zeroline = false;
    p.layout.xaxis1.showline = false;
    p.layout.yaxis1.showline = false;
    p.layout.xaxis1.showticklabels = false;
    p.layout.yaxis1.showticklabels = false;
    p.layout.xaxis1.domain = [0, 1];
    p.layout.yaxis1.domain = [0, 1];
    p.layout.plot_bgcolor = 'white';
    p.layout.paper_bgcolor = 'white';
    p.data{1}.colorbar.x = 1.02;
    p.data{1}.colorbar.xanchor = 'left';
    p.data{1}.colorbar.y = 0.5;
    p.data{1}.colorbar.yanchor = 'middle';
    p.data{1}.colorbar.len = 0.95;
    p.data{1}.colorbar.thickness = 20; % pixels
    p.data{1}.colorbar.thicknessmode = 'pixels';
    p.layout.margin.l = 0;
    p.layout.margin.r = 60;  % Adjust based on colorbar position
    p.layout.margin.t = 0;
    p.layout.margin.b = 0;
    p.layout.margin.pad = 0; 
    p.data{1}.z = round(p.data{1}.z,2);
    pause(3);
    progressMsg('cleanup_montages: completed...');
end


function [p] = cleanup_seg(p)
    progressMsg('cleanup_seg: starting...');

    p.layout.plot_bgcolor = 'black';
    p.layout.paper_bgcolor = 'black';
    p.data{1}.colorbar.tickcolor = "rgb(255,255,255)";
    p.data{1}.colorbar.tickfont.color = "rgb(255,255,255)";
    p.data{1}.colorbar.ticklabelposition = "outside";
    p.data{1}.colorbar.ticks = "outside";
    p.layout.autosize = true;
    for dd = 1 : 4
        p.data{dd}.z = round(p.data{dd}.z,2);
    end
    pause(3);
    progressMsg('cleanup_seg: completed...');
end


function [p] = cleanup_spectra(p,limits)
    progressMsg('cleanup_spectra: starting with %d traces...', length(p.data));

    for dd = 1 : length(p.data)
        progressMsg('cleanup_spectra: processing trace %d/%d, type = %s', dd, length(p.data), string(p.data{dd}.type));

        if strcmp(p.data{dd}.type,"heatmap")
            progressMsg('cleanup_spectra: rounding heatmap z data for trace %d...', dd);
            p.data{dd}.z = round(p.data{dd}.z,1);
        end

        if strcmp(p.data{dd}.type,"scatter") && length(p.data{dd}.x) == 4
           progressMsg('cleanup_spectra: converting 4-point scatter trace %d to line box...', dd);

           p.data{dd}.marker.size = 5;
           p.data{dd}.marker.line.width = 0;
           p.data{dd}.mode = 'lines';
           p.data{dd}.line.color = p.data{dd}.marker.color;
           p.data{dd}.line.width = 2.5;
           p.data{dd}.line.dash  = "solid";
           tempCorner3_x = p.data{dd}.x(3);
           tempCorner3_y = p.data{dd}.y(3);
           p.data{dd}.x(3)=p.data{dd}.x(4);
           p.data{dd}.y(3)=p.data{dd}.y(4);
           p.data{dd}.x(4)=tempCorner3_x;
           p.data{dd}.y(4)=tempCorner3_y;
           p.data{dd}.x(5)=p.data{dd}.x(1);
           p.data{dd}.y(5)=p.data{dd}.y(1);
        end

        if strcmp(p.data{dd}.type,"scatter") && length(p.data{dd}.x) > 50
            progressMsg('cleanup_spectra: rounding/cropping scatter trace %d with %d points...', dd, length(p.data{dd}.x));

            p.data{dd}.x = round(p.data{dd}.x,3);
            p.data{dd}.y = round(p.data{dd}.y,3);
            p.data{dd}.y = p.data{dd}.y(p.data{dd}.x > limits(1) & p.data{dd}.x < limits(2));
            p.data{dd}.x = p.data{dd}.x(p.data{dd}.x > limits(1) & p.data{dd}.x < limits(2));
        end
    end 
    pause(3);
    progressMsg('cleanup_spectra: completed...');
end


function progressMsg(varargin)
    fprintf(varargin{:});
    fprintf('\n');
    drawnow limitrate nocallbacks;
end

function [varStruct] = createInteractiveFigure(MRSCont,Module,varStruct,limits,outputFigures)
    switch Module
        case 'CheckPlotly'
            % Ensure that plotly works correctly
            progressMsg('Checking Plotly offline functionality...');
            try
              progressMsg('Determining report voxel indices...');
              if ~isfield(MRSCont.opts.MRSI.report,'VoxelIndices')
                  VoxelIndices = [round(MRSCont.raw{kk}.nXvoxels/2),round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                  round(MRSCont.raw{kk}.nXvoxels/2)+1,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                  round(MRSCont.raw{kk}.nXvoxels/2)+2,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);];
              else
                  VoxelIndices = MRSCont.opts.MRSI.report.VoxelIndices;
              end
            
              progressMsg('Plotly check: plotting raw spectra/localization...');
              osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyLoad',1,'Fit1DStack',0,2,1,0,0);
              fig = gcf;
              set(fig, 'Visible', 'off');
              set(fig,'Renderer','painters');
              drawnow limitrate nocallbacks; 
              PosLoadSpec = get(fig,'Position');
              set(gcf,'Position',[PosLoadSpec(1) PosLoadSpec(2) 4*PosLoadSpec(4) PosLoadSpec(4)])
            
              progressMsg('Plotly check: calling fig2plotly for Raw...');
              p = fig2plotly(gcf, 'offline', true,'filename','Raw','fileopt','new','open',false);
              close(fig);
              progressMsg('Plotly check completed successfully.');
            catch
              progressMsg('Plotly check failed. Installing plotly offline...');
              fprintf('Installing plotly for offline plotting.');
              plotlysetup_offline('http://cdn.plot.ly/plotly-latest.min.js');
              try
                progressMsg('Retrying Plotly export after installation...');
            
                if ~isfield(MRSCont.opts.MRSI.report,'VoxelIndices')
                    VoxelIndices = [round(MRSCont.raw{kk}.nXvoxels/2),round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                    round(MRSCont.raw{kk}.nXvoxels/2)+1,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                    round(MRSCont.raw{kk}.nXvoxels/2)+2,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);];
                else
                    VoxelIndices = MRSCont.opts.MRSI.report.VoxelIndices;
                end
            
                progressMsg('Fallback: plotting raw spectra/localization...');
                osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyLoad',1,'Fit1DStack',0,2,1,0,0);
                fig = gcf;
                set(fig, 'Visible', 'off');
                set(fig,'Renderer','painters');
                drawnow limitrate nocallbacks; 
                PosLoadSpec = get(fig,'Position');
                set(gcf,'Position',[PosLoadSpec(1) PosLoadSpec(2) 4*PosLoadSpec(4) PosLoadSpec(4)])
            
                progressMsg('Fallback: exporting Raw with fig2plotly...');
                p = fig2plotly(gcf, 'offline', true,'filename','Raw','fileopt','new','open',false);
            
                progressMsg('Fallback: cleaning Raw spectra...');
                p = cleanup_spectra(p,limits);
            
                progressMsg('Update plotly: calling plotly(p)...');
                plotly(p);
                progressMsg('Update plotly: completed.');
            
                progressMsg('Fallback: moving Raw.html...');
                movefile(fullfile(pwd,'Raw.html'),fullfile(outputFigures, 'Raw.html'));
                progressMsg('Moving file completed...');
                close(fig)
            
                progressMsg('Fallback Plotly export completed.');
              catch
                progressMsg('Fallback Plotly export failed.');
                fprintf('Failed to install plotly for offline HTML plotting. Consult getplotlyoffline.');
              end
            end
            progressMsg('Initial Plotly check finished. Closed all figures.');

        case 'OspreyLoadSpectra'
           
            % Plot spectra
            progressMsg('OspreyLoad: determining voxel indices...');
            if ~isfield(MRSCont.opts.MRSI.report,'VoxelIndices')
                VoxelIndices = [round(MRSCont.raw{kk}.nXvoxels/2),round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                round(MRSCont.raw{kk}.nXvoxels/2)+1,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                round(MRSCont.raw{kk}.nXvoxels/2)+2,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);];
            else
                VoxelIndices = MRSCont.opts.MRSI.report.VoxelIndices;
            end
           
            progressMsg('OspreyLoad: plotting raw spectra/localization...');
            osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyLoad',1,'Fit1DStack',0,2,1,0,0);
            
            fig = gcf;
            set(fig, 'Visible', 'off');
            set(fig,'Renderer','painters');
            drawnow limitrate nocallbacks;  
            
            varStruct.PosLoadSpec = get(fig,'Position');
            set(gcf,'Position',[varStruct.PosLoadSpec(1) varStruct.PosLoadSpec(2) 4*varStruct.PosLoadSpec(4) varStruct.PosLoadSpec(4)])
            
            progressMsg('OspreyLoad: exporting Raw.html with fig2plotly...');
            p = fig2plotly(fig, 'offline', true,'filename','Raw','fileopt','new','open',false);
        
            progressMsg('OspreyLoad: cleaning Raw spectra...');
            p = cleanup_spectra(p,limits);
        
            progressMsg('Update plotly: calling plotly(p)...');
            plotly(p);
            progressMsg('Update plotly: completed.');
            
            progressMsg('OspreyLoad: moving Raw.html to output folder...');
            movefile(fullfile(pwd,'Raw.html'),fullfile(outputFigures, 'Raw.html'));
            progressMsg('Moving file completed...');
            close(fig)
            progressMsg('OspreyLoad: Raw spectra export completed.');

        case 'OspreyLoadQuickMaps'
            % Quick Maps Load
            progressMsg('OspreyLoad: preparing quick maps...');
            currentFolder = pwd;
            field_names = fieldnames(MRSCont.opts.MRSI.quickMapsList{1});
            for ff = 1 : length(field_names)
                MRSCont.opts.MRSI.quickMaps.(field_names{ff}) = MRSCont.opts.MRSI.quickMapsList{1}.(field_names{ff});
            end
            
            progressMsg('OspreyLoad: starting raw quick map exports...');
            for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
                for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}))
                    spec = MRSCont.opts.MRSI.quickMaps.specs{ss};
                    name = MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll};
        
                    progressMsg('OspreyLoad: starting raw quick map %d/%d, %d/%d: %s %s', ...
                        ss, length(MRSCont.opts.MRSI.quickMaps.specs), ...
                        ll, length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})), ...
                        spec, name);
        
                    varStruct.names_quickMaps_raw{end+1} = [spec ' ' name];
                    osp_plotQuickmaps(MRSCont, spec, name);
        
                    progressMsg('OspreyLoad: completed raw quick map %d/%d, %d/%d: %s %s', ...
                        ss, length(MRSCont.opts.MRSI.quickMaps.specs), ...
                        ll, length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})), ...
                        spec, name);
        
                    fig = gcf;
                    set(fig, 'Visible', 'off');
                    set(fig,'Renderer','painters');
                    drawnow limitrate nocallbacks;
        
                    progressMsg('OspreyLoad: exporting quick map %s_%s.html...', spec, name);
                    p = fig2plotly(fig, 'offline', true,'filename',[spec '_' name],'fileopt','new','open',false);
        
                    progressMsg('OspreyLoad: cleaning quick map %s_%s...', spec, name);
                    p = cleanup_montages(p);
        
                    progressMsg('Update plotly: calling plotly(p)...');
                    plotly(p);
                    progressMsg('Update plotly: completed.');
        
                    progressMsg('OspreyLoad: moving quick map %s_%s.html...', spec, name);
                    movefile(fullfile(pwd,[spec '_' name '.html']),fullfile(outputFigures, [spec '_' name '.html']));
                    progressMsg('Moving file completed...');
                    varStruct.files_quickMaps_raw{end+1} = fullfile(outputFigures, [spec '_' name '.html']);
                    close(fig)
                    progressMsg('OspreyLoad: completed quick map %s_%s.html', spec, name);
                end
            end

        case 'OspreyCoreg'
                progressMsg('OspreyCoreg: plotting coregistration...');
                osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI', 1, 0, 0, 0,1, MRSCont.raw{1, 1}.nZvoxels);
                
                fig = gcf;
                set(fig, 'Visible', 'off');
                set(fig,'Renderer','painters');
                
                drawnow limitrate nocallbacks;  
                progressMsg('OspreyCoreg: exporting Coreg.png...');
                saveas(fig,'Coreg.png');
                
                varStruct.CoregPos = get(fig,'Position');
                
                progressMsg('OspreyCoreg: moving Coreg.png...');
                movefile(fullfile(pwd,'Coreg.png'),fullfile(outputFigures, 'Coreg.png'));
                progressMsg('Moving file completed...');
                close(fig)

        case 'OspreySeg'
            progressMsg('OspreySeg: plotting coregistration/segmentation overlay...');
            osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI', 0, 1, 1, 1,1, MRSCont.raw{1, 1}.nZvoxels);
            
            fig = gcf;
            set(fig, 'Visible', 'off');
            set(fig,'Renderer','painters');
            
            drawnow limitrate nocallbacks;  
            progressMsg('OspreySeg: exporting CoregSeg.png...');
            saveas(fig,'CoregSeg.png');
            
            varStruct.CoregSegPos = get(fig,'Position');
            
            progressMsg('OspreySeg: moving CoregSeg.png...');
            movefile(fullfile(pwd,'CoregSeg.png'),fullfile(outputFigures, 'CoregSeg.png'));
            progressMsg('Moving file completed...');
            close(fig)
            
            progressMsg('OspreySeg: plotting segmentation maps...');
            osp_plotSegMRSI(MRSCont,1,1,MRSCont.raw{1, 1}.nZvoxels, 1);
            
            fig = gcf;
            set(fig, 'Visible', 'off');
            set(fig,'Renderer','painters');
            
            drawnow limitrate nocallbacks;  
            set(fig, 'Units', 'Normalized')
            
            varStruct.PosAtlas = get(gcf,'OuterPosition');
            set(fig, 'OuterPosition', [varStruct.PosAtlas(1), varStruct.PosAtlas(2), 1/5 + 0.05, 0.96])
            set(fig, 'Units', 'Pixels')
            varStruct.PosAtlas = get(fig,'Position');
            
            progressMsg('OspreySeg: exporting Seg.html with fig2plotly...');
            p = fig2plotly(fig, 'offline', true,'filename','Seg','fileopt','new','open',false);
        
            progressMsg('OspreySeg: cleaning Seg plotly object...');
            p = cleanup_seg(p);
        
            progressMsg('Update plotly: calling plotly(p)...');
            plotly(p);
            progressMsg('Update plotly: completed.');
            
            progressMsg('OspreySeg: moving Seg.html...');
            movefile(fullfile(pwd,'Seg.html'),fullfile(outputFigures, 'Seg.html'));
            progressMsg('Moving file completed...');
            close(fig)

        case 'OspreyProcessSpectra'
            if ~isfield(MRSCont.opts.MRSI.report,'VoxelIndices')
                    VoxelIndices = [round(MRSCont.raw{kk}.nXvoxels/2),round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                    round(MRSCont.raw{kk}.nXvoxels/2)+1,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                    round(MRSCont.raw{kk}.nXvoxels/2)+2,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);];
                else
                    VoxelIndices = MRSCont.opts.MRSI.report.VoxelIndices;
            end

             % Plot spectra
            progressMsg('OspreyProcess: determining TargetSpec...');
            if isfield(MRSCont.processed,'A')
                varStruct.TargetSpec = 'A';
            end
            if isfield(MRSCont.processed,'AFID')
                varStruct.TargetSpec = 'AFID';
            end
            if isfield(MRSCont.processed,'diff1')
                varStruct.TargetSpec = 'diff1';
            end
            
            progressMsg('OspreyProcess: plotting processed spectra, TargetSpec = %s...', varStruct.TargetSpec);
            osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyProcess',varStruct.TargetSpec,'Fit1DStack',0,2,1,0,0);
            
            fig = gcf;
            set(fig, 'Visible', 'off');
            set(fig,'Renderer','painters');
            
            drawnow limitrate nocallbacks;  
            
            varStruct.PosProcSpec = get(fig,'Position');
            set(gcf,'Position',[varStruct.PosProcSpec(1) varStruct.PosProcSpec(2) 4*varStruct.PosProcSpec(4) varStruct.PosProcSpec(4)])
            
            progressMsg('OspreyProcess: exporting Process.html with fig2plotly...');
            p = fig2plotly(fig, 'offline', true,'filename','Process','fileopt','new','open',false);
        
            progressMsg('OspreyProcess: cleaning Process spectra...');
            p = cleanup_spectra(p,limits);
        
            progressMsg('Update plotly: calling plotly(p)...');
            plotly(p);
            progressMsg('Update plotly: completed.');
            
            progressMsg('OspreyProcess: moving Process.html...');
            movefile(fullfile(pwd,'Process.html'),fullfile(outputFigures, 'Process.html'));
            progressMsg('Moving file completed...');
            close(fig)

        case 'OspreyProcessQuickMaps'
            progressMsg('OspreyProcess: starting processed quick map exports...');
            for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
                for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}))
                    spec = MRSCont.opts.MRSI.quickMaps.specs{ss};
                    name = MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll};
        
                    progressMsg('OspreyProcess: starting quick map %d/%d, %d/%d: %s %s', ...
                        ss, length(MRSCont.opts.MRSI.quickMaps.specs), ...
                        ll, length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})), ...
                        spec, name);
        
                    varStruct.names_quickMaps_proc{end+1} = [spec ' ' name];
                    osp_plotQuickmaps(MRSCont, spec, name);
        
                    progressMsg('OspreyProcess: completed quick map %d/%d, %d/%d: %s %s', ...
                        ss, length(MRSCont.opts.MRSI.quickMaps.specs), ...
                        ll, length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})), ...
                        spec, name);
        
                    fig = gcf;
                    set(fig, 'Visible', 'off');
                    set(fig,'Renderer','painters');
                    drawnow limitrate nocallbacks;
                    p = fig2plotly(fig, 'offline', true,'filename',[spec '_' name],'fileopt','new','open',false);
                    p = cleanup_montages(p);
                    progressMsg('Update plotly: calling plotly(p)...');
                    plotly(p);
                    progressMsg('Update plotly: completed.');
                    movefile(fullfile(pwd,[spec '_' name '.html']),fullfile(outputFigures, [spec '_' name '.html']));
                    progressMsg('Moving file completed...');
                    varStruct.files_quickMaps_proc{end+1} = fullfile(outputFigures, [spec '_' name '.html']);
                    close(fig)
                    progressMsg('OspreyProcess: completed quick map %s_%s.html', spec, name);
                end
            end

        case 'OspreyProcessSNRMap'
             progressMsg('OspreyProcess: plotting A SNR map...');
            varStruct.names_quickMaps_proc{end+1} = ['A SNR'];
            osp_plotQuickmaps(MRSCont, 'A', 'SNR');
            fig = gcf;
            set(fig, 'Visible', 'off');   
            drawnow limitrate nocallbacks;
            p = fig2plotly(fig, 'offline', true,'filename',['A_SNR'],'fileopt','new','open',false);
            p = cleanup_montages(p);
            progressMsg('Update plotly: calling plotly(p)...');
            plotly(p);
            progressMsg('Update plotly: completed.');
            movefile(fullfile(pwd,'A_SNR.html'),fullfile(outputFigures,  'A_SNR.html'));
            progressMsg('Moving file completed...');
            varStruct.files_quickMaps_proc{end+1} = fullfile(outputFigures, 'A_SNR.html');
            close(fig)
            progressMsg('OspreyProcess: completed A_SNR.html');
        case 'OspreyProcessFWHMMap'
            progressMsg('OspreyProcess: plotting A FWHM map...');
            varStruct.names_quickMaps_proc{end+1} = ['A FWHM'];
            osp_plotQuickmaps(MRSCont, 'A', 'FWHM');
            fig = gcf;
            set(fig, 'Visible', 'off');
            set(fig,'Renderer','painters');
            drawnow limitrate nocallbacks;
            p = fig2plotly(fig, 'offline', true,'filename',['A_FWHM'],'fileopt','new','open',false);
            p = cleanup_montages(p);
            progressMsg('Update plotly: calling plotly(p)...');    
            plotly(p);
            progressMsg('Update plotly: completed.');
            movefile(fullfile(pwd,'A_FWHM.html'),fullfile(outputFigures,  'A_FWHM.html'));
            progressMsg('Moving file completed...');
            varStruct.files_quickMaps_proc{end+1} = fullfile(outputFigures, 'A_FWHM.html');
            close(fig)
            progressMsg('OspreyProcess: completed A_FWHM.html');
        case 'OspreyProcessGlobalQCMap'
            progressMsg('OspreyProcess: plotting Global QC filtering map...');
            varStruct.names_quickMaps_proc{end+1} = ['Global QC filtering'];
            osp_plotMetabolitemaps(MRSCont,'GlobalQC','tNAA_Acetyl_only',1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
            fig = gcf;
            set(fig, 'Visible', 'off');
            set(fig,'Renderer','painters');
            drawnow limitrate nocallbacks;
            p = fig2plotly(fig, 'offline', true,'filename',['Global_QC'],'fileopt','new','open',false);
            p = cleanup_montages(p);
            progressMsg('Update plotly: calling plotly(p)...');
            plotly(p);
            progressMsg('Update plotly: completed.');
            movefile(fullfile(pwd,'Global_QC.html'),fullfile(outputFigures,  'Global_QC.html'));
            progressMsg('Moving file completed...');
            varStruct.files_quickMaps_proc{end+1} = fullfile(outputFigures, 'Global_QC.html');
            close(fig)
            progressMsg('OspreyProcess: completed Global_QC.html');
            progressMsg('Finished OspreyProcess section.');
        case 'OspreyFit'
            if ~isfield(MRSCont.opts.MRSI.report,'VoxelIndices')
                    VoxelIndices = [round(MRSCont.raw{kk}.nXvoxels/2),round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                    round(MRSCont.raw{kk}.nXvoxels/2)+1,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                                    round(MRSCont.raw{kk}.nXvoxels/2)+2,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);];
                else
                    VoxelIndices = MRSCont.opts.MRSI.report.VoxelIndices;
            end

            progressMsg('OspreyFit: plotting fits...');
            osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,1,0,0);
            
            fig = gcf;
            set(fig, 'Visible', 'off');
            set(fig,'Renderer','painters');
            drawnow limitrate nocallbacks;  
            
            varStruct.PosFitSpec = get(fig,'Position');
            set(fig,'Position',[varStruct.PosFitSpec(1) varStruct.PosFitSpec(2) 4*varStruct.PosFitSpec(4) varStruct.PosFitSpec(4)])
            
            progressMsg('OspreyFit: exporting Fit.html with fig2plotly...');
            p = fig2plotly(gcf, 'offline', true,'filename','Fit','fileopt','new','open',false);
        
            progressMsg('OspreyFit: cleaning Fit spectra...');
            p = cleanup_spectra(p,limits);
        
            progressMsg('Update plotly: calling plotly(p)...');
            plotly(p);
            progressMsg('Update plotly: completed.');
            
            progressMsg('OspreyFit: moving Fit.html...');
            movefile(fullfile(pwd,'Fit.html'),fullfile(outputFigures, 'Fit.html'));
            progressMsg('Moving file completed...');
            close(fig)
        case 'OspreyQuantifyMap'

            if ~isfield(MRSCont.opts.MRSI.report, 'quantifications')
                quantifcations = {'tCr'};
            else
                quantifcations = MRSCont.opts.MRSI.report.quantifications;
            end
        
            if ~isfield(MRSCont.opts.MRSI.report, 'metabolites')
                metabolites = {'tNAA'};
            else
                metabolites = MRSCont.opts.MRSI.report.metabolites;
            end
        
            
            % Loop over quantifications
            progressMsg('OspreyQuantify: exporting metabolite maps without QC...');
            for qq = 1 : length(quantifcations)
                for mm = 1 : length(metabolites)
                    progressMsg('OspreyQuantify: Starting no-QC map %d/%d, %d/%d: %s %s', ...
                        qq, length(quantifcations), mm, length(metabolites), quantifcations{qq}, metabolites{mm});
                    osp_plotMetabolitemaps(MRSCont,quantifcations{qq},metabolites{mm},1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
                    progressMsg('OspreyQuantify: Completed no-QC map %d/%d, %d/%d: %s %s', ...
                        qq, length(quantifcations), mm, length(metabolites), quantifcations{qq}, metabolites{mm});
                    fig = gcf;
                    set(fig, 'Visible', 'off');
                    set(fig,'Renderer','painters');
                    drawnow limitrate nocallbacks;
                    p = fig2plotly(fig, 'offline', true,'filename',[quantifcations{qq},'_',metabolites{mm}],'fileopt','new','open',false);
                    p = cleanup_montages(p);
                    progressMsg('Update plotly: calling plotly(p)...');
                    plotly(p);
                    progressMsg('Update plotly: completed.');
                    movefile(fullfile(pwd,[quantifcations{qq},'_',metabolites{mm} '.html']),fullfile(outputFigures,  [quantifcations{qq},'_',metabolites{mm} '.html']));
                    progressMsg('Moving file completed...');
                    varStruct.files_Quantification_Maps{end+1} = fullfile(outputFigures, [quantifcations{qq},'_',metabolites{mm} '.html']);
                    varStruct.names_Quantification_Maps{end+1} = [quantifcations{qq},' ',metabolites{mm}];
                    close(fig)
                    progressMsg('OspreyQuantify: completed %s_%s.html', quantifcations{qq}, metabolites{mm});
                end
            end
        case 'OspreyQuantifyMapQC'
            if ~isfield(MRSCont.opts.MRSI.report, 'quantifications')
                quantifcations = {'tCr'};
            else
                quantifcations = MRSCont.opts.MRSI.report.quantifications;
            end
        
            if ~isfield(MRSCont.opts.MRSI.report, 'metabolites')
                metabolites = {'tNAA'};
            else
                metabolites = MRSCont.opts.MRSI.report.metabolites;
            end

            
            % Loop over quantifications
            progressMsg('OspreyQuantify: exporting QC maps...');
            for qq = 1 : length(quantifcations)
                for mm = 1 : length(metabolites)
                    progressMsg('OspreyQuantify: Starting QC map %d/%d, %d/%d: %s_QC %s', ...
                        qq, length(quantifcations), mm, length(metabolites), quantifcations{qq}, metabolites{mm});
                    osp_plotMetabolitemaps(MRSCont,[quantifcations{qq} '_QC'],metabolites{mm},1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
                    progressMsg('OspreyQuantify: Completed QC map %d/%d, %d/%d: %s_QC %s', ...
                        qq, length(quantifcations), mm, length(metabolites), quantifcations{qq}, metabolites{mm});
                    fig = gcf;
                    set(fig, 'Visible', 'off');
                    set(fig,'Renderer','painters');
                    drawnow limitrate nocallbacks;
                    p = fig2plotly(fig, 'offline', true,'filename',[quantifcations{qq},'_QC','_',metabolites{mm}],'fileopt','new','open',false);
                    p = cleanup_montages(p);
                    progressMsg('Update plotly: calling plotly(p)...');
                    plotly(p);
                    progressMsg('Update plotly: completed.');
                    movefile(fullfile(pwd,[quantifcations{qq},'_QC','_',metabolites{mm} '.html']),fullfile(outputFigures,  [quantifcations{qq},'_QC','_',metabolites{mm} '.html']));
                    progressMsg('Moving file completed...');
                    varStruct.files_Quantification_Maps{end+1} = fullfile(outputFigures, [quantifcations{qq},'_QC','_',metabolites{mm} '.html']);
                    varStruct.names_Quantification_Maps{end+1} = [quantifcations{qq},' ',metabolites{mm}, 'QC filter'];
                    close(fig)
                    progressMsg('OspreyQuantify: completed %s_QC_%s.html', quantifcations{qq}, metabolites{mm});
                end
            end
        case 'OspreyQuantifyMapQCFilt'
            if ~isfield(MRSCont.opts.MRSI.report, 'quantifications')
                quantifcations = {'tCr'};
            else
                quantifcations = MRSCont.opts.MRSI.report.quantifications;
            end
        
            if ~isfield(MRSCont.opts.MRSI.report, 'metabolites')
                metabolites = {'tNAA'};
            else
                metabolites = MRSCont.opts.MRSI.report.metabolites;
            end

            
            % Loop over quantifications
            progressMsg('OspreyQuantify: exporting QC-filtered maps...');
            for qq = 1 : length(quantifcations)
                for mm = 1 : length(metabolites)
                    progressMsg('OspreyQuantify: Starting QC-filtered map %d/%d, %d/%d: %s_QCfilt %s', ...
                        qq, length(quantifcations), mm, length(metabolites), quantifcations{qq}, metabolites{mm});
                    osp_plotMetabolitemaps(MRSCont,[quantifcations{qq} '_QCfilt'],metabolites{mm},1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
                    progressMsg('OspreyQuantify: Completed QC-filtered map %d/%d, %d/%d: %s_QCfilt %s', ...
                        qq, length(quantifcations), mm, length(metabolites), quantifcations{qq}, metabolites{mm});
                    fig = gcf;
                    set(fig, 'Visible', 'off');
                    set(fig,'Renderer','painters');
                    drawnow limitrate nocallbacks;
                    p = fig2plotly(fig, 'offline', true,'filename',[quantifcations{qq},'_QCfilt','_',metabolites{mm}],'fileopt','new','open',false);
                    p = cleanup_montages(p);
                    progressMsg('Update plotly: calling plotly(p)...');
                    plotly(p);
                    progressMsg('Update plotly: completed.');
                    movefile(fullfile(pwd,[quantifcations{qq},'_QCfilt','_',metabolites{mm} '.html']),fullfile(outputFigures,  [quantifcations{qq},'_QCfilt','_',metabolites{mm} '.html']));
                    progressMsg('Moving file completed...');
                    varStruct.files_Quantification_Maps_QCfilt{end+1} = fullfile(outputFigures, [quantifcations{qq},'_QCfilt','_',metabolites{mm} '.html']);
                    varStruct.names_Quantification_Maps_QCfilt{end+1} = [quantifcations{qq},' ',metabolites{mm}, 'QC filtered'];
                    close(fig)
                    progressMsg('OspreyQuantify: completed %s_QCfilt_%s.html', quantifcations{qq}, metabolites{mm});
                end
            end
        case 'OspreyQuantifyMapCRLB'
            if ~isfield(MRSCont.opts.MRSI.report, 'quantifications')
                quantifcations = {'tCr'};
            else
                quantifcations = MRSCont.opts.MRSI.report.quantifications;
            end
        
            if ~isfield(MRSCont.opts.MRSI.report, 'metabolites')
                metabolites = {'tNAA'};
            else
                metabolites = MRSCont.opts.MRSI.report.metabolites;
            end

            
            for mm = 1 : length(metabolites)
                progressMsg('OspreyQuantify: Starting CRLB map %d/%d: %s', mm, length(metabolites), metabolites{mm});
                osp_plotMetabolitemaps(MRSCont,'CRLBs',metabolites{mm},1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
                progressMsg('OspreyQuantify: Completed CRLB map %d/%d: %s', mm, length(metabolites), metabolites{mm});
                fig = gcf;
                set(fig, 'Visible', 'off');
                set(fig,'Renderer','painters');
                drawnow limitrate nocallbacks;
                p = fig2plotly(fig, 'offline', true,'filename',['CRLBs','_',metabolites{mm}],'fileopt','new','open',false);
                p = cleanup_montages(p);
                progressMsg('Update plotly: calling plotly(p)...');
                plotly(p);
                progressMsg('Update plotly: completed.');
                movefile(fullfile(pwd,['CRLBs','_',metabolites{mm} '.html']),fullfile(outputFigures,  ['CRLBs','_',metabolites{mm} '.html']));
                progressMsg('Moving file completed...');
                varStruct.files_CRLB_Maps{end+1} = fullfile(outputFigures, ['CRLBs','_',metabolites{mm} '.html']);
                varStruct.names_CRLB_Maps{end+1} = ['CRLBs', ' ',  quantifcations{qq},' ',metabolites{mm}, 'QC filtered'];
                close(fig)
                progressMsg('OspreyQuantify: completed CRLBs_%s.html', metabolites{mm});
            end
        case 'OspreyOverviewGlobalConc'
            quantifcations = MRSCont.opts.MRSI.GlobalConc.quantities;
            metabolites = MRSCont.opts.MRSI.GlobalConc.metabolites;
            
            for qq = 1 : length(quantifcations)
                for mm = 1 : length(metabolites)
                    progressMsg('OspreyOverview: GlobalConc %d/%d, %d/%d: %s %s', ...
                        qq, length(quantifcations), mm, length(metabolites), quantifcations{qq}, metabolites{mm});
                    
                    osp_plotGlobalConcentration(MRSCont,metabolites{mm},quantifcations{qq});
                    
                    fig = gcf;
                    set(fig, 'Visible', 'off');
                    set(fig,'Renderer','painters');
                    drawnow limitrate nocallbacks;  
                    set(fig, 'Units', 'Normalized')
                    
                    Pos = get(fig,'OuterPosition');
                    set(fig, 'OuterPosition', [Pos(1), Pos(2),0.4, 0.4])
                    set(fig, 'Units', 'Pixels')
                    varStruct.PosGlobalConc = get(fig,'Position');
                    saveas(fig, ['GlobalConc_' quantifcations{qq},'_',metabolites{mm} '.png']);
                    
                    movefile(fullfile(pwd,['GlobalConc_' quantifcations{qq},'_',metabolites{mm} '.png']),fullfile(outputFigures,['GlobalConc_' quantifcations{qq},'_',metabolites{mm} '.png']));
                    progressMsg('Moving file completed...');
                    varStruct.files_GlobalConc{end+1} = fullfile(outputFigures, ['GlobalConc_' quantifcations{qq},'_',metabolites{mm} '.png']);
                    varStruct.names_GlobalConc{end+1} = ['Global Concentration ', quantifcations{qq},' ',metabolites{mm}];
                    
                    close(fig)
                    progressMsg('OspreyOverview: completed GlobalConc_%s_%s.png', quantifcations{qq}, metabolites{mm});
                end
            end
        case 'OspreyOverviewAtlas'
            quantifcations = MRSCont.opts.MRSI.atlas.quantities;
            regions = MRSCont.opts.MRSI.report.atlasregion;
            
            for qq = 1 : length(quantifcations)
                for rr = 1 : length(regions)
                    progressMsg('OspreyOverview: Atlas result %d/%d, region %d/%d: %s %s', ...
                        qq, length(quantifcations), rr, length(regions), quantifcations{qq}, regions{rr});
                    
                    osp_plotInteractiveAtlasAnalysis(MRSCont,1,quantifcations{qq},MRSCont.opts.MRSI.report.metabolites,'T1w_rMRSI',regions{rr});
                    
                    fig = gcf;
                    set(fig, 'Visible', 'off');
                    set(fig,'Renderer','painters');
                    drawnow limitrate nocallbacks;  
                    set(fig, 'Units', 'Normalized')
                    
                    Pos = get(fig,'OuterPosition');
                    set(fig, 'OuterPosition', [Pos(1), Pos(2),0.4, 0.4])
                    set(fig, 'Units', 'Pixels')
                    varStruct.PosAtlas = get(fig,'Position');
                    saveas(fig, ['AtlasResults_' quantifcations{qq},'_region_',num2str(rr), '.png']);
                    
                    movefile(fullfile(pwd,['AtlasResults_' quantifcations{qq},'_region_',num2str(rr), '.png']),fullfile(outputFigures,['AtlasResults_' quantifcations{qq},'_region_',num2str(rr),'.png']));
                    progressMsg('Moving file completed...');
                    varStruct.files_atlas{end+1} = fullfile(outputFigures, ['AtlasResults_' quantifcations{qq},'_region_',num2str(rr), '.png']);
                    varStruct.names_atlas{end+1} = ['Atlas Results ', quantifcations{qq},' ',regions{rr}];
                    
                    close(fig)
                    progressMsg('OspreyOverview: completed AtlasResults_%s_region_%d.png', quantifcations{qq}, rr);
                end
            end

    end
end