clear; clc; close all;

rootDir = '/media/data3/Joanne_SRT_pw/'; 
% rootDir = '/home/aclexp/pinwei/Joanne_SRT_fMRI/'; 

batchFile = fullfile(rootDir, 'conn_out', 'v2_no_param.mat');
QCFile = fullfile(rootDir, 'conn_out', 'v2_no_param_qc.csv');

%% Quality Assurance
% https://web.conn-toolbox.org/fmri-methods/denoising-pipeline#h.p_BaXJei3yiEQh

global CONN_x;

dataValidityScore = conn_qascores('DataValidity', [], []);
dataQualityScore = conn_qascores('DataQuality', [], [], L2_COVARS, {});
sensitivityVars = {'QC_ProportionValidScans','QC_MeanMotion', ...
                   'QC_MeanGSchange','QC_NORM_struct', ...
                   'QC_DOF','QC_PeakFC','QC_StdFC'};
dataSensitivityScore = conn_qascores('DataSensitivity', [], [], sensitivityVars, 'extreme');
QCScoreTable = array2table( ...
    [dataValidityScore, dataQualityScore, dataSensitivityScore], ...
    'VariableNames', ['Validity', 'Quality', 'Sensitivity'] ...
);
writetable(QCScoreTable, QCFile);

conn_batch( ...
    'filename', batchFile, ...
    'QA.foldername', 'QA_scripted', ...
    'QA.plots', { ...
       'QA_NORM structural', ...     % (1) : structural data + outline of MNI TPM template
       'QA_NORM functional', ...     % (2) : mean functional data + outline of MNI TPM template
       'QA_NORM rois', ...           % (3) : ROI data + outline of MNI TPM template  
       'QA_REG functional', ...      % (10): display mean functional data + structural data overlay
       'QA_COREG functional', ...    % (7) : display same single-slice (z=0) across multiple sessions/datasets
       'QA_TIME functional', ...     % (8) : display same single-slice (z=0) across all timepoints within each session
       'QA_TIMEART functional', ...  % (9) : display same single-slice (z=0) across all timepoints within each session together with ART timeseries (global signal changes and framewise displacement)
       'QA_DENOISE histogram', ...   % (11) : histogram of voxel-to-voxel correlation values (before and after denoising)
       'QA_DENOISE timeseries', ...  % (12) : BOLD signal traces before and after denoising
       'QA_DENOISE FC-QC', ...       % (13) : histogram of FC-QC associations; between-subject correlation between QC (Quality Control) and FC (Functional Connectivity) measures
       'QA_DENOISE scatterplot', ... % (14) : scatterplot of FC (Functional Connectivity r coeff) vs. distance (mm)
       'QA_COV'  ...                 % (31) : histogram display of second-level variables
    } ...
)

% 'QA_REG structural', ...      % (4) : structural data + outline of ROI
% 'QA_REG functional', ...      % (5) : mean functional data + outline of ROI
% 'QA_REG mni', ...             % (6) : reference MNI structural template + outline of ROI
% 'QA_SPM design', ...          % (21) : SPM review design matrix (from SPM.mat files only)
% 'QA_SPM contrasts', ...       % (22) : SPM review contrast specification (from SPM.mat files only)
% 'QA_SPM results', ...         % (23) : SPM review contrast effect-size (from SPM.mat files only)

% batch.QA.l2covariates = ;
% batch.QA.conditions = ;
% batch.QA.controlcovariates = ;