clear; clc; close all;

%% Setup paths 

IP = 23;

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
batchFile = fullfile(rootDir, 'conn_out', 'no_PM_260423.mat');

[~, fn, ~] = fileparts(batchFile);
logFile = fullfile(rootDir, 'log', [fn, '.log']);

[fd, ~, ~] = fileparts(logFile);
if ~isfolder(fd), mkdir(fd); end

%% Configuration

SID_LIST = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25];
N_SUBJS  = length(SID_LIST);
N_RUNS   = 8; 
TR       = 2.0; 
SPACE    = 'MNI152NLin2009cAsym';

FWHM     = 6;            % smoothing fwhm (mm)

POLY_ORD = 3;            % polynomial detrending order
BP_HZ    = [0.008 0.09]; % band-pass filter (Hz)
SIMULT   = false;        % simultaneous regression & band-pass
MOT24    = true;         % add quadratic motion parameters 
N_ACOMP  = 20;           % number of aCompCor components
N_AROMA  = Inf;          % number of ICA-AROMA components

ADD_PM   = false;        % add parametric modulation

% Conditions
COND_OF_INTEREST = [ ...
    "str_r12", "str_r34", "str_r56", "str_r78", ... % "structured"
    "swi_r34", "swi_r56" ... % "switch"
];
COND_NAMES = ["random", COND_OF_INTEREST, "incorrect"];
N_CONDS = numel(COND_NAMES);

% Covariates
COVAR_NAMES = {'realignment', 'fd', 'scrubbing', 'aCompCor', 'aroma'};
if ADD_PM, COVAR_NAMES = [COVAR_NAMES, {'slope'}]; end
N_COVARS = numel(COVAR_NAMES);

L2_COVARS = {'QC_MeanMotion', 'QC_InvalidScans'};

% Confounds -- can be 'Grey Matter', 'White Matter', 'CSF', any ROI name, 
% any covariate name, or 'Effect of *' where * represents any condition name.
CONFOUND_NAMES = {'realignment', 'scrubbing'};
if N_ACOMP > 0, CONFOUND_NAMES = [CONFOUND_NAMES, {'aCompCor'}]; end
if N_AROMA > 0, CONFOUND_NAMES = [CONFOUND_NAMES, {'aroma'}]; end
CONFOUND_NAMES = [CONFOUND_NAMES, append('Effect of ', COND_NAMES)];

% ROI files
ROI_NAMES = {'atlas', 'networks'}; 
ROI_FILES = {fullfile(connDir, 'rois', 'atlas.nii'), ...
             fullfile(connDir, 'rois', 'networks.nii')}; 

% Name of analysis 
ANALYSIS_NAME = 'gPPI';

% Between-subject factors
BS_EFFECT_NAMES = {'AllSubjects', 'ExcludeOutlierSubjects'};
BS_CONTRAST = [1 0];

% Contrasts
CONTRASTS = struct(); 
for i = 1:numel(COND_OF_INTEREST)
    cond = COND_OF_INTEREST(i);
    CONTRASTS(i).saveas = strcat(cond, '_ran');
    CONTRASTS(i).between_conditions_names = {cond, 'random'};
    CONTRASTS(i).between_conditions_contrast = [1 -1];
end

CONTRASTS(i+1).saveas = 'str_main';
CONTRASTS(i+1).between_conditions_names = {'str_r12', 'str_r34', 'str_r56', 'str_r78', 'random'};
CONTRASTS(i+1).between_conditions_contrast = [1 1 1 1 -4];

CONTRASTS(i+2).saveas = 'swi_main';
CONTRASTS(i+2).between_conditions_names = {'swi_r34', 'swi_r56', 'random'};
CONTRASTS(i+2).between_conditions_contrast = [1 1 -2];

%% Start logging

diary(logFile); 

t = datetime('now'); 
fprintf('\r\n============ %s ============\r\n', t);
fprintf('Diary Start!\r\n');

%% Setup + Denoising

if exist(batchFile, "file")
    fprintf('\r\n%s exists. Skip setup and denoising ...\r\n', batchFile);

else
    fprintf('\r\nBuilding initial batch structure (Setup + Denoising) ...\r\n');
    fprintf('\r\nSmoothing FWHM (mm)                : %d', FWHM);
    fprintf('\r\nPolynomial detrending order        : %d', POLY_ORD);
    fprintf('\r\nBand-pass filter (Hz)              : %s', num2str(BP_HZ));
    fprintf('\r\nSimultaneous regression & band-pass: %s', log2str(SIMULT));
    fprintf('\r\nAdd quadratic motion parameters    : %s', log2str(MOT24));
    fprintf('\r\nNumber of aCompCor components      : %d', N_ACOMP);
    fprintf('\r\nNumber of ICA-AROMA components     : %d', N_AROMA);
    fprintf('\r\nAdd parametric modulation          : %s', log2str(ADD_PM));
    fprintf('\r\n\r\n');

    batch = struct();
    batch.filename = batchFile;
    batch.parallel.N = 0; % run locally

    % Setup section -------------------------------------------------------
    % https://web.conn-toolbox.org/fmri-methods/preprocessing-pipeline#h.p_tRw5CGulGxI9

    batch.Setup.isnew                 = 1;
    batch.Setup.nsubjects             = N_SUBJS;
    batch.Setup.nsessions             = N_RUNS; 
    batch.Setup.RT                    = TR;
    batch.Setup.conditions.names      = COND_NAMES;
    batch.Setup.conditions.param      = zeros(1, N_CONDS);
    batch.Setup.covariates.names      = COVAR_NAMES;
    batch.Setup.subjects.effect_names = L2_COVARS;

    batch.Setup.rois.names = ROI_NAMES;
    for iroi = 1:numel(ROI_NAMES)
        batch.Setup.rois.files{iroi} = ROI_FILES{iroi}; 
    end
    batch.Setup.rois.multiplelabels = 1;

    meanMotionValues  = zeros(N_SUBJS, N_RUNS);
    invalidScanCounts = zeros(N_SUBJS, N_RUNS);
    
    for isub = 1:N_SUBJS
        sid = SID_LIST(isub); 
        subj = sprintf('sub-%03d', sid);
    
        batch.Setup.structurals{isub} = fullfile(bidsDir, subj, 'anat', ...
            sprintf('sub-%03d_space-%s_desc-preproc_T1w.nii.gz', sid, SPACE));
        batch.Setup.masks.Grey{isub}  = fullfile(bidsDir, subj, 'anat', ...
            sprintf('sub-%03d_space-%s_label-GM_probseg.nii.gz', sid, SPACE));
        batch.Setup.masks.White{isub} = fullfile(bidsDir, subj, 'anat', ...
            sprintf('sub-%03d_space-%s_label-WM_probseg.nii.gz', sid, SPACE));
        batch.Setup.masks.CSF{isub}   = fullfile(bidsDir, subj, 'anat', ...
            sprintf('sub-%03d_space-%s_label-CSF_probseg.nii.gz', sid, SPACE));
        
        for irun = 1:N_RUNS
            batch.Setup.functionals{isub}{irun} = fullfile(bidsDir, subj, 'func', ...
                sprintf('sub-%03d_task-srttprob_run-%02d_space-%s_desc-preproc_bold.nii.gz', sid, irun, SPACE));
            
            % Condition event files
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
                
                batch.Setup.conditions.onsets{icond}{isub}{irun}    = condData.onset;
                batch.Setup.conditions.durations{icond}{isub}{irun} = condData.duration;
                
                % Including a parametric modulator
                if strcmp(condName, 'structured') && any(contains(COVAR_NAMES, 'slope'))
                    batch.Setup.conditions.param(icond) = find(strcmp(COVAR_NAMES, 'slope')); 
                end
            end
    
            % Covariate timeseries
            for icov = 1:N_COVARS
                covarName = COVAR_NAMES{icov};
                covarFile = fullfile(bidsDir, subj, 'func', ...
                    'covariates', sprintf('%s_run-%02d.tsv', covarName, irun));
                covarData = readtable(covarFile, 'FileType', 'text', 'Delimiter', '\t'); 
    
                if ~isempty(covarData)
                    batch.Setup.covariates.files{icov}{isub}{irun} = covarFile;
                end
    
                % Save for l2 covariates
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
    
    batch.Setup.preprocessing.fwhm  = FWHM;
    batch.Setup.preprocessing.steps = {'functional_smooth'};
    
    batch.Setup.outputfiles = [1 1 0 0 0 0]; % create confound beta-maps & confound-corrected timeseries
    batch.Setup.done        = 1;
    batch.Setup.overwrite   = 'No';

    % Denoising step ------------------------------------------------------
    % https://web.conn-toolbox.org/fmri-methods/denoising-pipeline
    
    batch.Denoising.confounds.names = CONFOUND_NAMES;
    batch.Denoising.detrending      = POLY_ORD;
    batch.Denoising.filter          = BP_HZ; 

    if SIMULT; batch.Denoising.regbp = 2; end 
    
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
    
    batch.Denoising.done      = 1;
    batch.Denoising.overwrite = 'No';

    t0 = tic;
    conn_batch(batch);
    fprintf('\r\n--- Initial Setup + Denoising done (Elapsed time: %.1f min) ---\r\n', toc(t0)/60);
end

%% Quality Assurance
% https://web.conn-toolbox.org/fmri-methods/denoising-pipeline#h.p_BaXJei3yiEQh

load(batchFile, 'CONN_x');
global CONN_x;

fprintf('\r\nComputing QC scores ...\r\n');
try
    t0 = tic;
    s1 = conn_qascores('DataValidity',    [], []);
    s2 = conn_qascores('DataQuality',     [], [], L2_COVARS, {});
    s3 = conn_qascores('DataSensitivity', [], [], [], [], 'extreme');
    mean_qc = mean([s1, s2, s3], 'omitnan');

    fprintf('\r\n--- QC scores calculation completed (Elapsed time: %.1f sec) ---\r\n', toc(t0));
    fprintf('\r\nData Validity score   : %.4f', s1);
    fprintf('\r\nData Quality score    : %.4f', s2);
    fprintf('\r\nData Sensitivity score: %.4f', s3);
    fprintf('\r\nMean QC score         : %.4f', mean_qc);

catch err
    fprintf('\r\nFailed to compute QC scores:\r\n%s\r\n', err.message);
    end_diary();
    exit
end

%% First-level analysis 
% https://web.conn-toolbox.org/fmri-methods/connectivity-measures

if contains({CONN_x.Analyses.name}, ANALYSIS_NAME)
    fprintf('\r\nAnalysis "%s" exists. Not overwriting ...\r\n', ANALYSIS_NAME);
else
    try
        fprintf('\r\nConducting "%s" analysis ...\r\n\r\n', ANALYSIS_NAME);
        t0 = tic;
        conn_batch( ...
            'filename', batchFile, ...
            'Analysis.name', ANALYSIS_NAME, ...
            'Analysis.measure', 3, ...     % regression (bivariate)
            'Analysis.modulation', 1, ...  % gPPI interaction effect
            'Analysis.conditions', COND_LIST, ... 
            'Analysis.type', 3, ...        % ROI-to-ROI & Seed-to-Voxel
            'Analysis.done', 1, ...
            'Analysis.overwrite', 'No' ...
        );
        fprintf('\r\n--- Analysis done (Elapsed time: %.1f min) ---\r\n', toc(t0)/60);

    catch err
        fprintf('\r\nFailed to run first-level analysis:\r\n%s\r\n', err.message);
        end_diary();
        exit
    end
end

%% Second-level analysis 

for i = 1:numel(CONTRASTS)
    C = CONTRASTS(i);
    fd = fullfile(rootDir, 'conn_out', 'no_PM_260421', 'results', 'secondlevel', ANALYSIS_NAME, C.saveas);
    
    if exist(fd, "dir")
        fprintf('\r\nContrast %s may has been analyzed. Skipping it ...\r\n', C.saveas);
    else
        try
            fprintf('\r\nAnalyzing contrast %s ...\r\n\r\n', C.saveas);
            t0 = tic;
            conn_batch( ...
                'filename', batchFile, ...
                'Results.saveas', C.saveas, ...
                'Results.name', ANALYSIS_NAME, ...
                'Results.between_subjects.effect_names', BS_EFFECT_NAMES, ...
                'Results.between_subjects.contrast', BS_CONTRAST, ...
                'Results.between_conditions.effect_names', C.between_conditions_names, ...
                'Results.between_conditions.contrast', C.between_conditions_contrast, ...
                'Results.display', 0, ...
                'Results.done', 1 ...
            );
            fprintf('\r\n--- Done (Elapsed time: %.1f min) ---\r\n', toc(t0)/60);
    
        catch err
            fprintf('\r\nFailed to analyze contrast "%s":\r\n%s\r\n', C.saveas, err.message);
            continue;
        end
    end
end

%% End logging

function end_diary()
    fprintf('\r\n\r\nDiary End!');
    t = datetime('now');  
    fprintf('\r\n====== %s ======\r\n', t);
    
    diary off;
end

function str = log2str(a)
    if a, str = 'true'; else, str = 'false'; end
end
