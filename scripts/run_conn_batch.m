clear; clc; close all;

%% Configuration

SID_LIST = [1, 2, 4];
N_SUBJS = length(SID_LIST);
N_RUNS = 8; 

TR = 2.0; % repetition time (seconds)
SPACE = 'MNI152NLin2009cAsym';

COND_NAMES = {'random', 'structured', 'switch', 'incorrect'};
N_CONDS = numel(COND_NAMES);
COVAR_NAMES = {'motion', 'global', 'aCompCor', 'behav'}; 
N_COVARS = numel(COVAR_NAMES);

FWHM = 6; % smoothing fwhm (mm)
OVERWRITE = false;

rootDir   = '/media/data3/Joanne_SRT_pw/'; 
bidsDir   = fullfile(rootDir, 'data', 'fmriprep_mini');
batchFile = fullfile(rootDir, 'derivatives', 'CONN', 'proj_mini_2.mat');

roiNames = {'custom', 'atlas', 'networks'}; 
% customAtlasFile = '/usr/local/abin/FS.afni.MNI2009c_asym.nii.gz'; % the FS-based AFNI-distributed atlas
customAtlasFile = '/media/data3/Joanne_SRT_pw/data/meta/FS.afni.MNI2009c_asym.nii.gz'
defaultAtlasFile = '/home/aclexp/mytools/matlab/conn/rois/atlas.nii';
defaultNetworksFile = '/home/aclexp/mytools/matlab/conn/rois/networks.nii';

%% Prepares batch structure
% https://github.com/alfnie/conn/blob/master/conn_batch.m

% load batchFile
batch.filename = batchFile;

%% Setup & Preprocessing steps
% https://web.conn-toolbox.org/fmri-methods/preprocessing-pipeline#h.p_tRw5CGulGxI9

batch.Setup.isnew = 1;
batch.Setup.nsubjects = N_SUBJS;
batch.Setup.nsessions = N_RUNS; 
batch.Setup.RT = TR;
batch.Setup.conditions.names = COND_NAMES;
batch.Setup.conditions.param = zeros(1, N_CONDS);
batch.Setup.covariates.names = COVAR_NAMES;

batch.Setup.rois.names = roiNames;
batch.Setup.rois.files{1} = customAtlasFile; 
batch.Setup.rois.files{2} = defaultAtlasFile;
batch.Setup.rois.files{3} = defaultNetworksFile;
batch.Setup.rois.multiplelabels = 1;

for isub = 1:numel(SID_LIST)
    sid = SID_LIST(isub); 
    subj = sprintf('sub-%03d', sid);

    % paths to structural images
    batch.Setup.structurals{isub} = fullfile( ...
        bidsDir, subj, 'anat', sprintf('sub-%03d_space-%s_desc-preproc_T1w.nii.gz', sid, SPACE));

    % paths to GM/WM/CSF masks
    batch.Setup.masks.Grey{isub} = fullfile( ...
        bidsDir, subj, 'anat', sprintf('sub-%03d_space-%s_label-GM_probseg.nii.gz', sid, SPACE));
    batch.Setup.masks.White{isub} = fullfile( ...
        bidsDir, subj, 'anat', sprintf('sub-%03d_space-%s_label-WM_probseg.nii.gz', sid, SPACE));
    batch.Setup.masks.CSF{isub} = fullfile( ...
        bidsDir, subj, 'anat', sprintf('sub-%03d_space-%s_label-CSF_probseg.nii.gz', sid, SPACE));
    
    for irun = 1:N_RUNS

        % paths to BOLD images
        batch.Setup.functionals{isub}{irun} = fullfile( ...
            bidsDir, subj, 'func', sprintf('sub-%03d_task-srttprob_run-%02d_space-%s_desc-preproc_bold.nii.gz', sid, irun, SPACE));
        
        % condition event files
        for icond = 1:N_CONDS
            condName = COND_NAMES{icond};
            condFile = fullfile( ...
                bidsDir, subj, 'func', 'events', sprintf('%s_run-%02d.tsv', condName, irun));
            condData = readtable(condFile, 'FileType', 'text', 'Delimiter', '\t'); 
            batch.Setup.conditions.onsets{icond}{isub}{irun} = condData.onset;
            batch.Setup.conditions.durations{icond}{isub}{irun} = condData.duration;
            if strcmp(condName, 'structured')
                batch.Setup.conditions.param(icond) = find(strcmp(COVAR_NAMES, 'behav')); 
            end
        end

        % covariate timeseries
        for icov = 1:N_COVARS
            covarName = COVAR_NAMES{icov};
            batch.Setup.covariates.files{icov}{isub}{irun} = fullfile( ...
                bidsDir, subj, 'func', 'covariates', sprintf('%s_run-%02d.tsv', covarName, irun));
        end
    end
end

batch.Setup.preprocessing.fwhm = FWHM;
batch.Setup.preprocessing.steps = {'functional_smooth'};

batch.Setup.outputfiles = [0,1,0,0,0,0]; % creates denoised output files
batch.Setup.done = 1;
batch.Setup.overwrite = 'Yes';

conn_batch(batch)

%% Denoising step

batch.Denoising.filter = [0.01, 0.10]; % band-pass filter (Hz)
batch.Denoising.detrending = 2;
batch.Denoising.confounds.names = ;
batch.Denoising.done = 1;
batch.Denoising.overwrite = 'Yes';

%% First-level analysis step

% batch.Analysis.name = ;
batch.Analysis.regressors = ;
batch.Analysis.measure = ;
batch.Analysis.modulation = ;
batch.Analysis.condition = ;
batch.Analysis.sources

batch.Analysis.done = 1;
batch.Analysis.overwrite = 'Yes';

%% Running the batch

conn_batch(batch)

% launches conn gui to explore results
% conn
% conn('load', batchFile);
% conn gui_results

% save(batchFile, 'batch')