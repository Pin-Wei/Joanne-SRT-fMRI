clear; clc; close all;

%% Setup paths 

IP = 37;

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
% batchFile = fullfile(rootDir, 'conn_out', 'no_PM_260423.mat');
batchFile = fullfile(rootDir, 'conn_out', 'temp.mat');

[~, fn, ~] = fileparts(batchFile);
logFile = fullfile(rootDir, 'logs', [fn, '.log']);

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
CONFOUND_NAMES = ['Grey Matter', CONFOUND_NAMES];

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

logFID = fopen(logFile, 'a'); 

t = datetime('now', 'Format', 'MMMM d, yyyy h:mm a'); 
fprintf(logFID, '\r\n======================== %s ========================', t);
fprintf(logFID, '\r\nStart Logging!');
fprintf(logFID, '\r\n');

%% Setup + Denoising

if exist(batchFile, "file")
    fprintf(logFID, '\r\n%s exists. Skip setup and denoising ...', batchFile);
    fprintf(logFID, '\r\n');
else
    fprintf(logFID, '\r\nBuilding initial batch structure (Setup + Denoising) ...');
    fprintf(logFID, '\r\n');
    fprintf(logFID, '\r\nSmoothing FWHM (mm)                : %d', FWHM);
    fprintf(logFID, '\r\nPolynomial detrending order        : %d', POLY_ORD);
    fprintf(logFID, '\r\nBand-pass filter (Hz)              : %s', num2str(BP_HZ));
    fprintf(logFID, '\r\nSimultaneous regression & band-pass: %s', log2str(SIMULT));
    fprintf(logFID, '\r\nAdd quadratic motion parameters    : %s', log2str(MOT24));
    fprintf(logFID, '\r\nNumber of aCompCor components      : %d', N_ACOMP);
    fprintf(logFID, '\r\nNumber of ICA-AROMA components     : %d', N_AROMA);
    fprintf(logFID, '\r\nAdd parametric modulation          : %s', log2str(ADD_PM));
    fprintf(logFID, '\r\n');

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
    fprintf(logFID, '\r\n--- Finished initial setup + denoising (Elapsed time: %.1f min) ---', toc(t0)/60);
    fprintf(logFID, '\r\n');
end

%% Quality Assurance
% https://web.conn-toolbox.org/fmri-methods/denoising-pipeline#h.p_BaXJei3yiEQh

try
    load(batchFile, 'CONN_x');
    global CONN_x;
    
    fprintf(logFID, '\r\nComputing QC scores ...');
    fprintf(logFID, '\r\n');
    t0 = tic;
    s1 = conn_qascores('DataValidity',    [], []);
    s2 = conn_qascores('DataQuality',     [], [], L2_COVARS, {});
    s3 = conn_qascores('DataSensitivity', [], [], [], [], 'extreme');
    mean_qc = mean([s1, s2, s3], 'omitnan');
    
    fprintf(logFID, '\r\nData Validity score   : %.4f', s1);
    fprintf(logFID, '\r\nData Quality score    : %.4f', s2);
    fprintf(logFID, '\r\nData Sensitivity score: %.4f', s3);
    fprintf(logFID, '\r\nMean QC score         : %.4f', mean_qc);
    fprintf(logFID, '\r\n');
    fprintf(logFID, '\r\n--- Finished QC scores calculation (Elapsed time: %.1f sec) ---\r\n', toc(t0));
    fprintf(logFID, '\r\n');

catch err
    fprintf(logFID, '\r\n*** Failed to compute QC scores ***\r\n');
    fprintf(logFID, '\r\n%s', err.message);
    end_logging(logFID);
    exit
end

%% First-level analysis 
% https://web.conn-toolbox.org/fmri-methods/connectivity-measures

if contains({CONN_x.Analyses.name}, ANALYSIS_NAME)
    fprintf(logFID, '\r\nAnalysis "%s" exists. Not overwriting ...\r\n', ANALYSIS_NAME);
    fprintf(logFID, '\r\n');
else
    try
        fprintf(logFID, '\r\nConducting "%s" analysis ...\r\n', ANALYSIS_NAME);
        fprintf(logFID, '\r\n');
        t0 = tic;
        conn_batch( ...
            'filename', batchFile, ...
            'Analysis.name', ANALYSIS_NAME, ...
            'Analysis.measure', 3, ...     % regression (bivariate)
            'Analysis.modulation', 1, ...  % gPPI interaction effect
            'Analysis.conditions', COND_NAMES, ... 
            'Analysis.type', 3, ...        % ROI-to-ROI & Seed-to-Voxel
            'Analysis.done', 1, ...
            'Analysis.overwrite', 'No' ...
        );
        fprintf(logFID, '\r\n--- Done (Elapsed time: %.1f min) ---', toc(t0)/60);
        fprintf(logFID, '\r\n');

    catch err
        fprintf(logFID, '\r\n*** Failed to run first-level analysis ***\r\n');
        fprintf(logFID, '\r\n%s', err.message);
        end_logging(logFID);
        exit
    end
end

%% Second-level analysis 

for i = 1:numel(CONTRASTS)
    C = CONTRASTS(i);
    fd = fullfile(rootDir, 'conn_out', 'no_PM_260421', 'results', 'secondlevel', ANALYSIS_NAME, C.saveas);
    
    if exist(fd, "dir")
        fprintf(logFID, '\r\nContrast "%s" may has been analyzed. Skipping it ...', C.saveas);
        fprintf(logFID, '\r\n');
    else
        try
            fprintf(logFID, '\r\nAnalyzing contrast "%s" ...', C.saveas);
            fprintf(logFID, '\r\n');
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
            fprintf(logFID, '\r\n--- Done (Elapsed time: %.1f min) ---', toc(t0)/60);
            fprintf(logFID, '\r\n');

        catch err
            fprintf(logFID, '\r\n*** Failed to analyze contrast "%s" ***\r\n', C.saveas);
            fprintf(logFID, '\r\n%s', err.message);
            continue;
        end
    end
end

end_logging(logFID);

%% Functions

function end_logging(fid)
    fprintf(fid, '\r\n');
    fprintf(fid, '\r\nEnd Logging!');
    t = datetime('now', 'Format', 'MMMM d, yyyy h:mm a'); 
    fprintf(fid, '\r\n======================== %s ========================', t);
    fprintf(fid, '\r\n');
    fclose(fid);
end

function str = log2str(a)
    if a, str = 'true'; else, str = 'false'; end
end
