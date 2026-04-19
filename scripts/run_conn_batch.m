clear; clc; close all;

%% Configuration

IP = 37; % 23

if IP == 37
    rootDir = '/media/data3/Joanne_SRT_pw/';
    connDir = '/home/aclexp/mytools/matlab/conn';
elseif IP == 23
    rootDir = '/home/aclexp/pinwei/Joanne_SRT_fMRI/';
    connDir = '/home/aclexp/Software/conn';
else 
    error('Root directory for this IP address has not yet been defined.');
end

bidsDir = fullfile(rootDir, 'data', 'fmriprep');
batchFile = fullfile(rootDir, 'conn_out', 'v0_no_param.mat');

SID_LIST = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25];
N_SUBJS = length(SID_LIST);
N_RUNS = 8; 

TR = 2.0; % repetition time (seconds)
SPACE = 'MNI152NLin2009cAsym';

FWHM = 6; % smoothing fwhm (mm)
POLY_ORD = 3; % polynomial detrending order
BP_HZ = [0.008 0.09]; % band-pass filter (Hz)
SIMULT = true;  
MOT24 = false; 
N_ACOMP = 15;
N_AROMA = Inf; % 0

% COND_NAMES = {'random', 'structured', 'switch', 'incorrect'};
COND_NAMES = {'random', 'str_r12', 'str_r34', 'str_r56', 'str_r78', 'swi_r34', 'swi_r56', 'incorrect'};
N_CONDS = numel(COND_NAMES);
 
COVAR_NAMES = {'realignment', 'fd', 'scrubbing', 'aCompCor', 'aroma'}; 
% COVAR_NAMES = [COVAR_NAMES, {'slope'}];
N_COVARS = numel(COVAR_NAMES);

L2_COVARS = {'QC_MeanMotion', 'QC_InvalidScans'};

CONFOUND_NAMES = {'realignment', 'scrubbing'};
if N_ACOMP > 0, CONFOUND_NAMES = [CONFOUND_NAMES, {'aCompCor'}]; end
if N_AROMA > 0, CONFOUND_NAMES = [CONFOUND_NAMES, {'aroma'}]; end
CONFOUND_NAMES = [CONFOUND_NAMES, append('Effect of ', COND_NAMES)];
    % can be 'Grey Matter', 'White Matter', 'CSF', any ROI name, 
    % any covariate name, or 'Effect of *' where * represents any condition name.

ROI_NAMES = {'atlas', 'networks'}; 
ROI_FILES = {fullfile(connDir, 'rois', 'atlas.nii'), fullfile(connDir, 'rois', 'networks.nii')}; 

%% Prepares batch structure

batch.filename = batchFile;
batch.parallel.N = 0; % run locally

%% Setup & Preprocessing steps
% https://web.conn-toolbox.org/fmri-methods/preprocessing-pipeline#h.p_tRw5CGulGxI9

batch.Setup.isnew = 1;

batch.Setup.nsubjects = N_SUBJS;
batch.Setup.nsessions = N_RUNS; 
batch.Setup.RT = TR;
batch.Setup.conditions.names = COND_NAMES;
batch.Setup.conditions.param = zeros(1, N_CONDS);
batch.Setup.covariates.names = COVAR_NAMES;
batch.Setup.subjects.effect_names = L2_COVARS;

batch.Setup.rois.names = ROI_NAMES;
for iroi = 1:numel(ROI_NAMES)
    batch.Setup.rois.files{iroi} = ROI_FILES{iroi}; 
end
batch.Setup.rois.multiplelabels = 1;

meanMotionValues = zeros(N_SUBJS, N_RUNS);
invalidScanCounts = zeros(N_SUBJS, N_RUNS);

for isub = 1:N_SUBJS
    sid = SID_LIST(isub); 
    subj = sprintf('sub-%03d', sid);

    % paths to structural images
    batch.Setup.structurals{isub} = fullfile(bidsDir, subj, 'anat', ...
        sprintf('sub-%03d_space-%s_desc-preproc_T1w.nii.gz', sid, SPACE));

    % paths to GM/WM/CSF masks
    batch.Setup.masks.Grey{isub} = fullfile(bidsDir, subj, 'anat', ...
        sprintf('sub-%03d_space-%s_label-GM_probseg.nii.gz', sid, SPACE));
    batch.Setup.masks.White{isub} = fullfile(bidsDir, subj, 'anat', ...
        sprintf('sub-%03d_space-%s_label-WM_probseg.nii.gz', sid, SPACE));
    batch.Setup.masks.CSF{isub} = fullfile(bidsDir, subj, 'anat', ...
        sprintf('sub-%03d_space-%s_label-CSF_probseg.nii.gz', sid, SPACE));
    
    for irun = 1:N_RUNS

        % paths to BOLD images
        batch.Setup.functionals{isub}{irun} = fullfile(bidsDir, subj, 'func', ...
            sprintf('sub-%03d_task-srttprob_run-%02d_space-%s_desc-preproc_bold.nii.gz', sid, irun, SPACE));
        
        % condition event files
        for icond = 1:N_CONDS
            condName = COND_NAMES{icond};

            if contains(condName, 'str_')
                if ~contains(condName, num2str(irun))
                    batch.Setup.conditions.onsets{icond}{isub}{irun} = [];
                    batch.Setup.conditions.durations{icond}{isub}{irun} = [];
                    continue
                end
                condName = 'structured';
            elseif contains(condName, 'swi_')
                if ~contains(condName, num2str(irun))
                    batch.Setup.conditions.onsets{icond}{isub}{irun} = [];
                    batch.Setup.conditions.durations{icond}{isub}{irun} = [];
                    continue
                end
                condName = 'switch';
            end
            
            condFile = fullfile(bidsDir, subj, 'func', ...
                'events', sprintf('%s_run-%02d.tsv', condName, irun));
            condData = readtable(condFile, 'FileType', 'text', 'Delimiter', '\t');
            batch.Setup.conditions.onsets{icond}{isub}{irun} = condData.onset;
            batch.Setup.conditions.durations{icond}{isub}{irun} = condData.duration;
            
            % including a parametric modulator
            if strcmp(condName, 'structured') && any(contains(COVAR_NAMES, 'slope'))
                batch.Setup.conditions.param(icond) = find(strcmp(COVAR_NAMES, 'slope')); 
            end
        end

        % covariate timeseries
        for icov = 1:N_COVARS
            covarName = COVAR_NAMES{icov};
            covarFile = fullfile(bidsDir, subj, 'func', ...
                'covariates', sprintf('%s_run-%02d.tsv', covarName, irun));
            covarData = readtable(covarFile, 'FileType', 'text', 'Delimiter', '\t'); 
            if ~isempty(covarData)
                batch.Setup.covariates.files{icov}{isub}{irun} = covarFile;
            end

            % save for l2 covariates
            if matches(covarName, 'scrubbing')
                invalidScanCounts(isub, irun) = sum(covarData{:, 1});
            elseif matches(covarName, 'fd')
                meanMotionValues(isub, irun) = mean(covarData{:, 1});
            end
        end
    end
end

% l2 covariates (for QC)
for icov2 = 1:numel(L2_COVARS)
    covarName = L2_COVARS{icov2};
    if matches(covarName, 'QC_InvalidScans')
        batch.Setup.subjects.effects{icov2} = sum(invalidScanCounts, 2);
    elseif matches(covarName, 'QC_MeanMotion')
        batch.Setup.subjects.effects{icov2} = mean(meanMotionValues, 2);
    end
end

batch.Setup.preprocessing.fwhm = FWHM;
batch.Setup.preprocessing.steps = {'functional_smooth'};

batch.Setup.outputfiles = [1,1,0,0,0,0]; 
    % 1) creates confound beta-maps
    % 2) creates confound-corrected timeseries
    % 3) creates seed-to-voxel r-maps
    % 4) creates seed-to-voxel p-maps
    % 5) creates seed-to-voxel FDR-p-maps
    % 6) creates ROI-extraction REX files

batch.Setup.done = 1;
batch.Setup.overwrite = 'No';

%% Denoising step
% https://web.conn-toolbox.org/fmri-methods/denoising-pipeline

batch.Denoising.confounds.names = CONFOUND_NAMES;
batch.Denoising.detrending = POLY_ORD;
batch.Denoising.filter = BP_HZ; 
if SIMULT; batch.Denoising.regbp = 2; end % simultaneous regression & band-pass

idx = find(strcmp(CONFOUND_NAMES, 'realignment'));
batch.Denoising.confounds.deriv{idx} = 1;
if MOT24; batch.Denoising.confounds.power{idx} = 2; end

if any(contains(CONFOUND_NAMES, 'aCompCor'))
    idx = find(strcmp(CONFOUND_NAMES, 'aCompCor'));
    batch.Denoising.confounds.dimensions{idx} = N_ACOMP;
end

if any(contains(CONFOUND_NAMES, 'aroma'))
    idx = find(strcmp(CONFOUND_NAMES, 'aroma'));
    batch.Denoising.confounds.dimensions{idx} = N_AROMA;
end

for conf = ["White Matter", "CSF"]
    if any(contains(CONFOUND_NAMES, conf))
        idx = find(strcmp(CONFOUND_NAMES, conf));
        batch.Denoising.confounds.dimensions{idx} = 5;
    end
end

batch.Denoising.done = 1;
batch.Denoising.overwrite = 'No';

%% Finally

conn_batch(batch)

