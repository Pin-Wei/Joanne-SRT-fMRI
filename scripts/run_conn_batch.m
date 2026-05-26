function run_conn_batch(varargin)
    clc; close all;
    
    % Configuration -------------------------------------------------------
    
    p = inputParser;
    addRequired(p, 'IP', @(x) any(x == [23, 37])); % positional arguments
    addRequired(p, 'PROJECT', @ischar);
    parse(p, varargin{:});

    IP = p.Results.IP;
    PROJECT = p.Results.PROJECT;
    
    FWHM     = 6;            % smoothing fwhm (mm)
    POLY_ORD = 3;            % polynomial detrending order
    BP_HZ    = [0.008 0.09]; % band-pass filter (Hz)
    SIMULT   = false;        % simultaneous regression & band-pass
    MOT24    = true;         % add quadratic motion parameters 
    N_ACOMP  = 5;            % number of aCompCor components
    N_AROMA  = 10;           % number of ICA-AROMA components
    RM_GMR   = true;         % remove the average signal within the grey matter (GM) mask
    
    % Conditions and Contrasts
    [COND_OF_INTEREST, CONTRASTS] = define_for_version(PROJECT);

    COND_NAMES = ["random", COND_OF_INTEREST, "incorrect"];
    N_CONDS = numel(COND_NAMES);
    
    % Covariates
    COVAR_NAMES = {'realignment', 'fd', 'scrubbing', 'aCompCor', 'aroma'};
    N_COVARS = numel(COVAR_NAMES);

    L2_COVARS = {'QC_MeanMotion', 'QC_InvalidScans'};
    
    % Confounds
    CONFOUND_NAMES = {'realignment', 'scrubbing'};
    if N_ACOMP > 0, CONFOUND_NAMES = [CONFOUND_NAMES, {'aCompCor'}]; end
    if N_AROMA > 0, CONFOUND_NAMES = [CONFOUND_NAMES, {'aroma'}]; end
    CONFOUND_NAMES = [CONFOUND_NAMES, append('Effect of ', COND_NAMES)]; 
    if RM_GMR, CONFOUND_NAMES = ['Grey Matter', CONFOUND_NAMES]; end
    
    % Name of analysis 
    ANALYSIS_NAME = 'gPPI';
    
    % Between-subject factors
    BS_EFFECT_NAMES = {'AllSubjects', 'ExcludeOutlierSubjects'};
    BS_CONTRAST = [1 0];
    
    % Fixed information
    SID_LIST = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25];
    N_SUBJS  = length(SID_LIST);
    N_RUNS   = 8; 
    TR       = 2.0; 
    SPACE    = 'MNI152NLin2009cAsym';

    % Paths 
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
    batchFile = fullfile(rootDir, 'conn_out', [PROJECT, '.mat']);
    logFile = fullfile(rootDir, 'logs', [PROJECT, '.log']);
    
    [fd, ~, ~] = fileparts(logFile);
    if ~isfolder(fd), mkdir(fd); end
    
    % Atlases
    ROI_NAMES = {'atlas', 'networks', 'FS', 'BA'}; 
    ROI_FILES = { ...
        fullfile(connDir, 'rois', 'atlas.nii'), ...
        fullfile(connDir, 'rois', 'networks.nii'), ...
        fullfile(rootDir, 'data', 'meta', 'FS_atlas.nii'), ...
        fullfile(rootDir, 'data', 'meta', 'Brodmann_atlas.nii') ...
    }; 
    
    % Start logging -------------------------------------------------------
    
    logFID = fopen(logFile, 'a'); 
    t = datetime('now', 'Format', 'MMMM d, yyyy h:mm a'); 
    write_log(logFID, '\r\n======================== %s ========================', t);
    write_log(logFID, '\r\nStart Logging!');
    write_log(logFID, '\r\n');
        
    if exist(batchFile, 'file')
        write_log(logFID, '\r\n"%s" exists.', batchFile);
        write_log(logFID, '\r\nSkip setup and denoising ...');
        write_log(logFID, '\r\n');
    else
        write_log(logFID, '\r\nBuilding initial batch structure (Setup + Denoising) ...');
        write_log(logFID, '\r\n');
        write_log(logFID, '\r\nSmoothing FWHM (mm)                : %d', FWHM);
        write_log(logFID, '\r\nPolynomial detrending order        : %d', POLY_ORD);
        write_log(logFID, '\r\nBand-pass filter (Hz)              : %s', num2str(BP_HZ));
        write_log(logFID, '\r\nSimultaneous regression & band-pass: %s', log2str(SIMULT));
        write_log(logFID, '\r\nAdd quadratic motion parameters    : %s', log2str(MOT24));
        write_log(logFID, '\r\nNumber of aCompCor components      : %d', N_ACOMP);
        write_log(logFID, '\r\nNumber of ICA-AROMA components     : %d', N_AROMA);
        write_log(logFID, '\r\n');
    
        % Setup -----------------------------------------------------------
        % https://web.conn-toolbox.org/fmri-methods/preprocessing-pipeline#h.p_tRw5CGulGxI9

        batch = struct();
        batch.filename = batchFile;
        batch.parallel.N = 0; % run locally

        batch.Setup.isnew                 = 1;
        batch.Setup.nsubjects             = N_SUBJS;
        batch.Setup.nsessions             = N_RUNS; 
        batch.Setup.RT                    = TR;
        batch.Setup.conditions.names      = COND_NAMES;
        batch.Setup.conditions.param      = zeros(1, N_CONDS);
        batch.Setup.covariates.names      = COVAR_NAMES;
        batch.Setup.subjects.effect_names = L2_COVARS;
        batch.Setup.subjects.group_names  = {'AllSubjects'};
        batch.Setup.subjects.groups       = ones(N_SUBJS, 1);
        batch.Setup.subjects.add          = 1;
    
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
                    [onsets, durations] = cond_setter(condName, irun, subj, bidsDir);
                    batch.Setup.conditions.onsets{icond}{isub}{irun} = onsets;
                    batch.Setup.conditions.durations{icond}{isub}{irun} = durations;
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
        batch.Setup.overwrite   = 1;
    
        % Denoising -------------------------------------------------------
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
        batch.Denoising.overwrite = 1;
    
        t0 = tic;
        conn_batch(batch);
        write_log(logFID, '\r\n--- Finished initial setup + denoising (Elapsed time: %s) ---', format_time(toc(t0)));
        write_log(logFID, '\r\n');
    end
    
    % Quality Assurance ---------------------------------------------------
    % https://web.conn-toolbox.org/fmri-methods/denoising-pipeline#h.p_BaXJei3yiEQh
    
    try
        conn('load', batchFile); 
        global CONN_x;
    
        [fp, fn, ~] = fileparts(batchFile);
        globQA = dir(fullfile(fp, fn, 'results', 'qa', 'QA_GUIrequest_*'));
    
        if ~isempty(globQA)
            write_log(logFID, '\r\nQC scores may have been computed. Skipping ...');
            write_log(logFID, '\r\n');
        else
            write_log(logFID, '\r\nComputing QC scores ...');
            write_log(logFID, '\r\n');
            t0 = tic;
            s1 = conn_qascores('DataValidity', [], []);
            s2 = conn_qascores('DataQuality', [], [], L2_COVARS, {});
            s3 = conn_qascores('DataSensitivity', [], [], [], [], 'extreme');
            mean_qc = mean([s1, s2, s3], 'omitnan');
            conn save; 
    
            write_log(logFID, '\r\nData Validity score   : %.4f', s1);
            write_log(logFID, '\r\nData Quality score    : %.4f', s2);
            write_log(logFID, '\r\nData Sensitivity score: %.4f', s3);
            write_log(logFID, '\r\nMean QC score         : %.4f', mean_qc);
            write_log(logFID, '\r\n');
            write_log(logFID, '\r\n--- Finished QC scores calculation (Elapsed time: %s) ---', format_time(toc(t0)));
            write_log(logFID, '\r\n');
        end
    catch err
        write_log(logFID, '\r\n*** Failed to compute QC scores ***');
        write_log(logFID, '\r\n%s', err.message);
        write_log(logFID, '\r\n');
        end_logging(logFID);
        exit
    end
    
    % First-level analysis ------------------------------------------------
    % https://web.conn-toolbox.org/fmri-methods/connectivity-measures
    
    if any(matches({CONN_x.Analyses.name}, ANALYSIS_NAME))
        write_log(logFID, '\r\nAnalysis "%s" exists. Not overwriting ...', ANALYSIS_NAME);
        write_log(logFID, '\r\n');
    else
        try
            write_log(logFID, '\r\nConducting "%s" analysis ...', ANALYSIS_NAME);
            write_log(logFID, '\r\n');
            t0 = tic;
            conn_batch( ...
                'filename', batchFile, ...
                'Analysis.name', ANALYSIS_NAME, ...
                'Analysis.measure', 3, ...        % regression (bivariate)
                'Analysis.modulation', 1, ...     % gPPI interaction effect
                'Analysis.conditions', COND_NAMES, ... 
                'Analysis.type', 3, ...           % ROI-to-ROI & Seed-to-Voxel
                ... % 'Analysis.sources', ANA_SOURCES, ...
                'Analysis.done', 1, ...
                'Analysis.overwrite', 1 ...
            );
            conn save; 
            write_log(logFID, '\r\n--- Done (Elapsed time: %s) ---', format_time(toc(t0)));
            write_log(logFID, '\r\n');
    
        catch err
            write_log(logFID, '\r\n*** Failed to run first-level analysis ***');
            write_log(logFID, '\r\n%s', err.message);
            write_log(logFID, '\r\n');
            end_logging(logFID);
            exit
        end
    end
    
    % Second-level analysis -----------------------------------------------
    
    for i = 1:numel(CONTRASTS)
        C = CONTRASTS(i);
        
        if isfolder(fullfile(CONN_x.folders.secondlevel, ANALYSIS_NAME, C.saveas))
        % if any(matches(CONN_x.Results.saved.names, C.saveas))
            write_log(logFID, '\r\nContrast "%s" may has been analyzed. Skipping ...', C.saveas);
            write_log(logFID, '\r\n');
        else
            try
                write_log(logFID, '\r\nAnalyzing contrast "%s" ...', C.saveas);
                write_log(logFID, '\r\n');
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
                conn save; 
                write_log(logFID, '\r\n--- Done (Elapsed time: %s) ---', format_time(toc(t0)));
                write_log(logFID, '\r\n');
    
            catch err
                write_log(logFID, '\r\n*** Failed to analyze contrast "%s" ***', C.saveas);
                write_log(logFID, '\r\n%s', err.message);
                write_log(logFID, '\r\n');
                continue;
            end
        end
    end
    
    end_logging(logFID);
end

%% Other functions

function [COND_OF_INTEREST, CONTRASTS] = define_for_version(PROJECT)
    if ismember(PROJECT, {'no_PM_260504', 'no_PM_260515'})
        COND_OF_INTEREST = [ ...
            "str_r12", "str_r34", "str_r56", "str_r78", ... % "structured"
            "swi_r34", "swi_r56" ... % "switch"
        ];
        CONTRASTS = init_contrast(COND_OF_INTEREST);
        CONTRASTS = add_contrast(CONTRASTS, ...
            'str_main', ...
            {'str_r12', 'str_r34', 'str_r56', 'str_r78', 'random'}, ...
            [1 1 1 1 -4] ...
        );
        CONTRASTS = add_contrast(CONTRASTS, ...
            'swi_main', ...
            {'swi_r34', 'swi_r56', 'random'}, ...
            [1 1 -2] ...
        );

    elseif ismember(PROJECT, {'no_PM_260506', 'no_PM_260512', 'no_PM_260524'})
        COND_OF_INTEREST = [ ...
            "str_fst_r1", "str_fst_r2", "str_fst_r3", "str_fst_r4", "str_fst_r5", "str_fst_r6", ... 
            "str_snd_r3", "str_snd_r4", "str_snd_r5", "str_snd_r6", "str_snd_r7", "str_snd_r8", ... 
            "swi_r3", "swi_r4", "swi_r5", "swi_r6" ... 
        ];
        CONTRASTS = init_contrast(COND_OF_INTEREST);
        CONTRASTS = add_contrast(CONTRASTS, ...
            'str_fst', ...
            {'str_fst_r1', 'str_fst_r2', 'str_fst_r3', ...
             'str_fst_r4', 'str_fst_r5', 'str_fst_r6', 'random'}, ...
            [1 1 1 1 1 1 -6] ...
        );
        CONTRASTS = add_contrast(CONTRASTS, ...
            'str_snd', ...
            {'str_snd_r3', 'str_snd_r4', 'str_snd_r5', ...
             'str_snd_r6', 'str_snd_r7', 'str_snd_r8', 'random'}, ...
            [1 1 1 1 1 1 -6] ...
        ); 
        CONTRASTS = add_contrast(CONTRASTS, ...
            'switch', ...
            {'swi_r3', 'swi_r4', 'swi_r5', 'swi_r6', 'random'}, ...
            [1 1 1 1 -4] ...
        ); 
    else
        error("Error! Project not defined.")
    end
end

function CONTRASTS = init_contrast(COND_OF_INTEREST)
    CONTRASTS = struct(); 
    for i = 1:numel(COND_OF_INTEREST)
        cond = COND_OF_INTEREST(i);
        CONTRASTS(i).saveas = strcat(cond, '_ran');
        CONTRASTS(i).between_conditions_names = {cond, 'random'};
        CONTRASTS(i).between_conditions_contrast = [1 -1];
    end
end

function CONTRASTS = add_contrast(CONTRASTS, name, conds, weights)
    CONTRASTS(end+1).saveas = name;
    CONTRASTS(end).between_conditions_names = conds;
    CONTRASTS(end).between_conditions_contrast = weights;
end

function [onsets, durations] = cond_setter(condName, irun, subj, bidsDir)
    onsets = []; durations = []; 

    prefixes = {'str_', 'str_fst', 'str_snd', 'swi_'};
    if any(startsWith(condName, prefixes)) && ~contains(condName, num2str(irun))
        return; 
    end

    if startsWith(condName, 'str_fst')
        fileTag = 'str_fst';
    elseif startsWith(condName, 'str_snd')
        fileTag = 'str_snd';
    elseif startsWith(condName, 'str_')
        fileTag = 'structured';
    elseif startsWith(condName, 'swi_')
        fileTag = 'switch';
    else
        fileTag = condName; 
    end
    condFile  = fullfile(bidsDir, subj, 'func', ...
        'events', sprintf('%s_run-%02d.tsv', fileTag, irun));

    if isfile(condFile)
        condData  = readtable(condFile, 'FileType', 'text', 'Delimiter', '\t');
        onsets    = condData.onset; 
        durations = condData.duration;
    else
        warning('File not found: %s', condFile);
    end
end

function write_log(fid, fmt, varargin)
    msg  = sprintf(fmt, varargin{:});
    fprintf('%s', msg);
    if fid > 0, fprintf(fid, '%s', msg); end
end

function end_logging(fid)
    write_log(fid, '\r\nEnd Logging!');
    t = datetime('now', 'Format', 'MMMM d, yyyy h:mm a'); 
    write_log(fid, '\r\n======================== %s ========================', t);
    write_log(fid, '\r\n');
    fclose(fid);
end

function str = log2str(a)
    if a, str = 'true'; 
    else, str = 'false'; 
    end
end

function str = format_time(t)
    if t/3600 > 0
        hr = floor(t / 3600);
        min = floor(mod(t, 3600) / 60);
        sec = mod(t, 60);
        str = sprintf('%d hr %d min %.1f sec', hr, min, sec);
    elseif t/60 > 0
        min = floor(t / 60);
        sec = mod(t, 60);
        str = sprintf('%d min %.1f sec', min, sec);
    else
        str = sprintf('%.1f sec', t);
    end
end