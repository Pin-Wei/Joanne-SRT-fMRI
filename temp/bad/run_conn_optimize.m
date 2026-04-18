%% run_conn_optimize.m
% =========================================================================
% Two-stage grid-search optimization of CONN fMRI denoising parameters.
% Finds the parameter combination that maximizes CONN QC scores
% (Data Validity, Data Quality, Data Sensitivity).
%
% Approach
% --------
%   Stage 1 (coarse) : tests a reduced grid of representative levels
%                      (3 x 2 x 2 x 2 x 3 x 2 = 144 combinations)
%   Stage 2 (fine)   : fixes the winning discrete parameters and sweeps
%                      all N_ACOMP and N_AROMA levels around the winner
%
% Speed trick: after one-time Setup + Denoising (to extract confound
% timeseries into COV files), subsequent iterations only modify
% CONN_x.Preproc in memory and call conn_qascores, which re-applies
% denoising on the sample voxels on-the-fly.  This is orders of
% magnitude faster than re-running full voxel-level denoising.
%
% Outputs (all in tempDir)
% -------
%   optimization_results.xlsx  - every combination tested, sorted by
%                                Mean_QC descending; second sheet has
%                                only the best pipeline
%   optimization_report.txt    - human-readable summary with best combo,
%                                main-effect analysis, recommendation
%   optimization_log.txt       - timestamped progress log for monitoring
%                                (tail -f optimization_log.txt)
%
% Prerequisites
% -------------
%   - CONN toolbox + SPM on the MATLAB path
%   - fMRIPrep outputs at bidsDir with pre-extracted covariate TSVs
%
% Usage
%   1. Edit the Configuration section below (paths, SKIP_SETUP flag)
%   2. >> run_conn_optimize
%   3. Monitor progress: tail -f .../temp/optimization_log.txt
%
% NOTE: The original run_conn_preproc.m is never modified.
% =========================================================================

clear; clc; close all;
global CONN_x;

%% ====================== 1. CONFIGURATION ================================

rootDir  = '/home/aclexp/pinwei/Joanne_SRT_fMRI/';
bidsDir  = fullfile(rootDir, 'data', 'fmriprep');
tempDir  = fullfile(rootDir, 'temp');
optBatch = fullfile(rootDir, 'conn_out', 'optimize_denoise.mat');

logFile    = fullfile(tempDir, 'optimization_log.txt');
xlsFile    = fullfile(tempDir, 'optimization_results.xlsx');
reportFile = fullfile(tempDir, 'optimization_report.txt');

% --- Fixed experimental parameters (not optimized) -----------------------
SID_LIST = [1,2,4,6,8,9,10,12,13,14,15,16,17,18,19,20,21,22,24,25];
N_SUBJS  = numel(SID_LIST);
N_RUNS   = 8;
TR       = 2.0;
SPACE    = 'MNI152NLin2009cAsym';
FWHM     = 6;

COND_NAMES  = {'random','str_r12','str_r34','str_r56','str_r78', ...
               'swi_r34','swi_r56','incorrect'};
COVAR_NAMES = {'realignment','fd','scrubbing','aCompCor','aroma'};
L2_COVARS   = {'QC_MeanMotion','QC_InvalidScans'};
SENS_VARS   = {'QC_ProportionValidScans','QC_MeanMotion', ...
               'QC_MeanGSchange','QC_NORM_struct', ...
               'QC_DOF','QC_PeakFC','QC_StdFC'};

ROI_NAMES = {'atlas','networks'};
ROI_FILES = {'/home/aclexp/Software/conn/rois/atlas.nii', ...
             '/home/aclexp/Software/conn/rois/networks.nii'};

% --- Full parameter search space ------------------------------------------
POLY_VALS   = [1, 2, 3];                       % polynomial detrending
BP_OPTS     = {[0.008, 0.09]; [0.008, Inf]};   % band-pass / high-pass
SIMULT_VALS = [1, 0];                          % simultaneous regression
MOT24_VALS  = [1, 0];                          % 24-param motion model
ACOMP_ALL   = [5, 10, 15, 20, 30, 40, 50];    % aCompCor components
AROMA_ALL   = [0, 5, 10, Inf];                 % ICA-AROMA components

% --- Coarse-grid subsets for Stage 1 --------------------------------------
ACOMP_COARSE = [5, 20, 50];
AROMA_COARSE = [0, Inf];

% --- Control flags --------------------------------------------------------
% Set SKIP_SETUP = true if the optimize_denoise project already exists
% from a previous run (skips the slow Setup+Denoising init)
SKIP_SETUP = false;

%% ====================== 2. INITIAL SETUP ================================

logfid = fopen(logFile, 'w');
log_msg(logfid, '====== CONN Denoising Optimization Started ======');
log_msg(logfid, 'MATLAB %s  |  CONN %s', version, conn('ver'));
log_msg(logfid, 'Project: %s', optBatch);

if ~SKIP_SETUP
    log_msg(logfid, 'Building initial batch structure (Setup + Denoising)...');

    batch = struct();
    batch.filename   = optBatch;
    batch.parallel.N = 0;

    % --- Setup section (mirrors run_conn_preproc.m) -----------------------
    batch.Setup.isnew      = 1;
    batch.Setup.nsubjects  = N_SUBJS;
    batch.Setup.nsessions  = N_RUNS;
    batch.Setup.RT         = TR;
    batch.Setup.conditions.names = COND_NAMES;
    batch.Setup.conditions.param = zeros(1, numel(COND_NAMES));
    batch.Setup.covariates.names = COVAR_NAMES;
    batch.Setup.subjects.effect_names = L2_COVARS;

    batch.Setup.rois.names = ROI_NAMES;
    for iroi = 1:numel(ROI_NAMES)
        batch.Setup.rois.files{iroi} = ROI_FILES{iroi};
    end
    batch.Setup.rois.multiplelabels = 1;

    meanMotionValues  = zeros(N_SUBJS, N_RUNS);
    invalidScanCounts = zeros(N_SUBJS, N_RUNS);

    for isub = 1:N_SUBJS
        sid  = SID_LIST(isub);
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
                sprintf('sub-%03d_task-srttprob_run-%02d_space-%s_desc-preproc_bold.nii.gz', ...
                sid, irun, SPACE));

            for icond = 1:numel(COND_NAMES)
                condName = COND_NAMES{icond};
                if contains(condName, 'str_')
                    if ~contains(condName, num2str(irun))
                        batch.Setup.conditions.onsets{icond}{isub}{irun}    = [];
                        batch.Setup.conditions.durations{icond}{isub}{irun} = [];
                        continue
                    end
                    condName = 'structured';
                elseif contains(condName, 'swi_')
                    if ~contains(condName, num2str(irun))
                        batch.Setup.conditions.onsets{icond}{isub}{irun}    = [];
                        batch.Setup.conditions.durations{icond}{isub}{irun} = [];
                        continue
                    end
                    condName = 'switch';
                end
                condFile = fullfile(bidsDir, subj, 'func', 'events', ...
                    sprintf('%s_run-%02d.tsv', condName, irun));
                condData = readtable(condFile, 'FileType','text', 'Delimiter','\t');
                batch.Setup.conditions.onsets{icond}{isub}{irun}    = condData.onset;
                batch.Setup.conditions.durations{icond}{isub}{irun} = condData.duration;
            end

            for icov = 1:numel(COVAR_NAMES)
                covarName = COVAR_NAMES{icov};
                covarFile = fullfile(bidsDir, subj, 'func', 'covariates', ...
                    sprintf('%s_run-%02d.tsv', covarName, irun));
                covarData = readtable(covarFile, 'FileType','text', 'Delimiter','\t');
                if ~isempty(covarData)
                    batch.Setup.covariates.files{icov}{isub}{irun} = covarFile;
                end
                if matches(covarName, 'scrubbing')
                    invalidScanCounts(isub, irun) = sum(covarData{:,1});
                elseif matches(covarName, 'fd')
                    meanMotionValues(isub, irun)  = mean(covarData{:,1});
                end
            end
        end
    end

    for icov2 = 1:numel(L2_COVARS)
        if matches(L2_COVARS{icov2}, 'QC_InvalidScans')
            batch.Setup.subjects.effects{icov2} = sum(invalidScanCounts, 2);
        elseif matches(L2_COVARS{icov2}, 'QC_MeanMotion')
            batch.Setup.subjects.effects{icov2} = mean(meanMotionValues, 2);
        end
    end

    batch.Setup.preprocessing.fwhm  = FWHM;
    batch.Setup.preprocessing.steps = {'functional_smooth'};
    batch.Setup.outputfiles = [1,1,0,0,0,0];
    batch.Setup.done      = 1;
    batch.Setup.overwrite  = 'No';

    % --- Initial Denoising (uses original defaults) -----------------------
    INIT_CONFOUNDS = [{'realignment','scrubbing','aCompCor','aroma'}, ...
                      strcat('Effect of ', COND_NAMES)];
    batch.Denoising.confounds.names = INIT_CONFOUNDS;
    batch.Denoising.detrending = 1;
    batch.Denoising.filter     = [0.008 Inf];
    batch.Denoising.regbp      = 2;  % simultaneous

    idx = find(strcmp(INIT_CONFOUNDS, 'realignment'));
    batch.Denoising.confounds.deriv{idx} = 1;
    batch.Denoising.confounds.power{idx} = 2;
    idx = find(strcmp(INIT_CONFOUNDS, 'aCompCor'));
    batch.Denoising.confounds.dimensions{idx} = max(ACOMP_ALL);  % extract maximum
    idx = find(strcmp(INIT_CONFOUNDS, 'aroma'));
    batch.Denoising.confounds.dimensions{idx} = Inf;             % all components
    batch.Denoising.done      = 1;
    batch.Denoising.overwrite = 'Yes';

    log_msg(logfid, 'Running conn_batch (Setup + initial Denoising) ...');
    t_setup = tic;
    conn_batch(batch);
    log_msg(logfid, 'Initial Setup + Denoising done (%.1f min)', toc(t_setup)/60);
else
    log_msg(logfid, 'SKIP_SETUP=true: loading existing project ...');
    conn('load', optBatch);
    log_msg(logfid, 'Project loaded');
end

% Verify CONN_x is populated
assert(~isempty(CONN_x) && isfield(CONN_x,'Preproc'), ...
    'CONN_x not initialised.  Set SKIP_SETUP=false and re-run.');

% Save the fully-initialised Preproc as our template
orig_preproc = CONN_x.Preproc;
log_msg(logfid, 'Saved orig_preproc template (%d confounds)', ...
    numel(orig_preproc.confounds.names));

%% ====================== 3. OPTIMISATION LOOP ============================
%
% Pre-allocate result storage.
%   Columns: [iter, stage, poly, bp_idx, bp_lo, bp_hi,
%             simult, mot24, acomp, aroma,
%             validity, quality, sensitivity, mean_qc, runtime_s]
%
MAX_ITER = 500;
R = nan(MAX_ITER, 15);
S = cell(MAX_ITER, 1);   % status string per iteration
iter = 0;                 % running counter

% ---- Stage 1: coarse grid -----------------------------------------------
log_msg(logfid, '');
log_msg(logfid, '========== STAGE 1: COARSE GRID ==========');

[g1,g2,g3,g4,g5,g6] = ndgrid(1:numel(POLY_VALS),   1:numel(BP_OPTS), ...
                               1:numel(SIMULT_VALS),  1:numel(MOT24_VALS), ...
                               1:numel(ACOMP_COARSE), 1:numel(AROMA_COARSE));
grid = [g1(:), g2(:), g3(:), g4(:), g5(:), g6(:)];
n_combos = size(grid, 1);
log_msg(logfid, 'Stage 1: %d combinations', n_combos);

for ic = 1:n_combos
    poly   = POLY_VALS(grid(ic,1));
    bp     = BP_OPTS{grid(ic,2)};       bp_idx = grid(ic,2);
    simult = SIMULT_VALS(grid(ic,3));
    mot24  = MOT24_VALS(grid(ic,4));
    acomp  = ACOMP_COARSE(grid(ic,5));
    aroma  = AROMA_COARSE(grid(ic,6));

    iter = iter + 1;
    log_msg(logfid, '[S1 %03d/%03d] POLY=%d  BP=[%.3f,%.3g]  SIM=%d  M24=%d  AC=%d  AR=%g', ...
        ic, n_combos, poly, bp(1), bp(2), simult, mot24, acomp, aroma);

    [s1, s2, s3, mean_qc, elapsed, status] = run_one_combo( ...
        orig_preproc, poly, bp, simult, mot24, acomp, aroma, ...
        L2_COVARS, SENS_VARS);

    R(iter,:) = [iter, 1, poly, bp_idx, bp(1), bp(2), ...
                 simult, mot24, acomp, aroma, s1, s2, s3, mean_qc, elapsed];
    S{iter}   = status;
    log_msg(logfid, '  -> V=%.4f  Q=%.4f  S=%.4f  Mean=%.4f  (%.1fs)  [%s]', ...
        s1, s2, s3, mean_qc, elapsed, status);
end

% ---- Identify Stage 1 winner --------------------------------------------
s1_rows  = find([R(1:iter,2)] == 1);
[~, best_s1_local] = max(R(s1_rows, 14));   % column 14 = mean_qc
best_s1  = s1_rows(best_s1_local);
win_poly   = R(best_s1, 3);
win_bp_idx = R(best_s1, 4);
win_simult = R(best_s1, 7);
win_mot24  = R(best_s1, 8);
win_acomp  = R(best_s1, 9);
win_aroma  = R(best_s1, 10);

log_msg(logfid, '');
log_msg(logfid, 'Stage 1 WINNER (iter %d): POLY=%d  BP_idx=%d  SIM=%d  M24=%d  AC=%d  AR=%g  Mean_QC=%.4f', ...
    best_s1, win_poly, win_bp_idx, win_simult, win_mot24, win_acomp, win_aroma, R(best_s1,14));

% ---- Stage 2: fine grid around winner ------------------------------------
log_msg(logfid, '');
log_msg(logfid, '========== STAGE 2: FINE GRID ==========');

% Determine N_ACOMP neighbourhood: all values within +/-15 of winner
acomp_fine = ACOMP_ALL(abs(ACOMP_ALL - win_acomp) <= 15);
% Always include the immediate neighbours from the full list
win_pos = find(ACOMP_ALL == win_acomp, 1);
lo = max(1, win_pos - 1);  hi = min(numel(ACOMP_ALL), win_pos + 1);
acomp_fine = unique([acomp_fine, ACOMP_ALL(lo:hi)]);

% All AROMA levels
aroma_fine = AROMA_ALL;

% Build Stage 2 grid (fix poly/bp/simult/mot24, sweep acomp x aroma)
[ga, gr] = ndgrid(1:numel(acomp_fine), 1:numel(aroma_fine));
grid2 = [ga(:), gr(:)];

% Remove combinations already tested in Stage 1
already = false(size(grid2,1), 1);
for j = 1:size(grid2,1)
    ac = acomp_fine(grid2(j,1));
    ar = aroma_fine(grid2(j,2));
    already(j) = any(R(s1_rows,9)==ac & R(s1_rows,10)==ar & ...
                     R(s1_rows,3)==win_poly & R(s1_rows,4)==win_bp_idx & ...
                     R(s1_rows,7)==win_simult & R(s1_rows,8)==win_mot24);
end
grid2 = grid2(~already, :);
n_fine = size(grid2, 1);
log_msg(logfid, 'Stage 2: %d new combinations (acomp_fine=%s, aroma_fine=%s)', ...
    n_fine, mat2str(acomp_fine), mat2str(aroma_fine));

bp_win = BP_OPTS{win_bp_idx};

for ic = 1:n_fine
    acomp = acomp_fine(grid2(ic,1));
    aroma = aroma_fine(grid2(ic,2));

    iter = iter + 1;
    log_msg(logfid, '[S2 %03d/%03d] POLY=%d  BP=[%.3f,%.3g]  SIM=%d  M24=%d  AC=%d  AR=%g', ...
        ic, n_fine, win_poly, bp_win(1), bp_win(2), win_simult, win_mot24, acomp, aroma);

    [s1, s2, s3, mean_qc, elapsed, status] = run_one_combo( ...
        orig_preproc, win_poly, bp_win, win_simult, win_mot24, acomp, aroma, ...
        L2_COVARS, SENS_VARS);

    R(iter,:) = [iter, 2, win_poly, win_bp_idx, bp_win(1), bp_win(2), ...
                 win_simult, win_mot24, acomp, aroma, s1, s2, s3, mean_qc, elapsed];
    S{iter}   = status;
    log_msg(logfid, '  -> V=%.4f  Q=%.4f  S=%.4f  Mean=%.4f  (%.1fs)  [%s]', ...
        s1, s2, s3, mean_qc, elapsed, status);
end

% Trim to actual number of iterations
R = R(1:iter, :);
S = S(1:iter);

log_msg(logfid, '');
log_msg(logfid, 'All %d iterations complete.', iter);

%% ====================== 4. EXPORT RESULTS ===============================

log_msg(logfid, 'Exporting results ...');

% ---- 4a. Build results table ---------------------------------------------
% Map BP index and AROMA to human-readable strings
BP_str  = cell(iter, 1);
AR_str  = cell(iter, 1);
for k = 1:iter
    BP_str{k} = sprintf('[%.3f, %s]', R(k,5), num2str(R(k,6)));
    if isinf(R(k,10)), AR_str{k} = 'Inf'; else, AR_str{k} = num2str(R(k,10)); end
end

T = table( ...
    (1:iter)',         R(:,2),          R(:,3), ...
    BP_str,            logical(R(:,7)), logical(R(:,8)), ...
    R(:,9),            R(:,10), ...
    R(:,11),           R(:,12),         R(:,13),          R(:,14), ...
    S,                 R(:,15), ...
    'VariableNames', {'Iter','Stage','POLY_ORD', ...
                      'BP_HZ','SIMULT','MOT24', ...
                      'N_ACOMP','N_AROMA', ...
                      'Validity','Quality','Sensitivity','Mean_QC', ...
                      'Status','Runtime_s'});

% Sort by Mean_QC descending (NaN at end)
[~, sort_idx] = sort(T.Mean_QC, 'descend', 'MissingPlacement','last');
T = T(sort_idx, :);

% ---- 4b. Write XLSX Sheet 1: all results --------------------------------
writetable(T, xlsFile, 'Sheet','All Results');

% ---- 4c. Write XLSX Sheet 2: best pipeline -------------------------------
best_row = T(1, :);
bp_best  = BP_OPTS{best_row.Stage};  % note: use actual bp from best row
summary_text = sprintf([ ...
    'Best denoising pipeline found by grid search optimization.\n' ...
    'POLY_ORD = %d  |  BP_HZ = %s  |  SIMULT = %d  |  MOT24 = %d\n' ...
    'N_ACOMP = %d  |  N_AROMA = %s\n' ...
    'Validity = %.4f  |  Quality = %.4f  |  Sensitivity = %.4f  |  Mean_QC = %.4f\n' ...
    'Tested %d parameter combinations across 2 stages.'], ...
    best_row.POLY_ORD, best_row.BP_HZ{1}, best_row.SIMULT, best_row.MOT24, ...
    best_row.N_ACOMP, AR_str{sort_idx(1)}, ...
    best_row.Validity, best_row.Quality, best_row.Sensitivity, best_row.Mean_QC, ...
    iter);

best_table = best_row;
writetable(best_table, xlsFile, 'Sheet','Best Pipeline');
% Append summary as a second small table
writetable(cell2table({summary_text}, 'VariableNames', {'Summary'}), ...
    xlsFile, 'Sheet','Best Pipeline', 'Range','A4');

log_msg(logfid, 'Results written to %s', xlsFile);

% ---- 4d. Write summary report -------------------------------------------
rfid = fopen(reportFile, 'w');

fprintf(rfid, '# Denoising Pipeline Optimization Report\n');
fprintf(rfid, '# Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(rfid, '# Project:   %s\n', optBatch);
fprintf(rfid, '# Subjects:  %d  |  Runs: %d\n\n', N_SUBJS, N_RUNS);

% ---- Best combination ---------------------------------------------------
fprintf(rfid, '## Best Parameter Combination\n\n');
fprintf(rfid, '  %-16s  %s\n', 'Parameter', 'Value');
fprintf(rfid, '  %-16s  %s\n', '----------------', '----------');
fprintf(rfid, '  %-16s  %d\n',   'POLY_ORD',  best_row.POLY_ORD);
fprintf(rfid, '  %-16s  %s\n',   'BP_HZ',     best_row.BP_HZ{1});
fprintf(rfid, '  %-16s  %d\n',   'SIMULT',    best_row.SIMULT);
fprintf(rfid, '  %-16s  %d\n',   'MOT24',     best_row.MOT24);
fprintf(rfid, '  %-16s  %d\n',   'N_ACOMP',   best_row.N_ACOMP);
fprintf(rfid, '  %-16s  %s\n\n', 'N_AROMA',   AR_str{sort_idx(1)});
fprintf(rfid, '  %-16s  %.4f\n', 'Validity',    best_row.Validity);
fprintf(rfid, '  %-16s  %.4f\n', 'Quality',     best_row.Quality);
fprintf(rfid, '  %-16s  %.4f\n', 'Sensitivity', best_row.Sensitivity);
fprintf(rfid, '  %-16s  %.4f\n\n','Mean_QC',    best_row.Mean_QC);

% ---- Score distributions (main effects) ----------------------------------
fprintf(rfid, '## Score Distributions (main-effect summary)\n\n');
fprintf(rfid, 'For each parameter, the mean QC score across all tested\n');
fprintf(rfid, 'combinations that used that level (marginal average).\n\n');

valid_mask = ~isnan(R(:,14));  % only successful runs
R_valid = R(valid_mask, :);

param_info = { ...
    'POLY_ORD',  3,  POLY_VALS;
    'BP_HZ',     4,  1:numel(BP_OPTS);
    'SIMULT',    7,  SIMULT_VALS;
    'MOT24',     8,  MOT24_VALS;
    'N_ACOMP',   9,  unique(R_valid(:,9))';
    'N_AROMA',  10,  unique(R_valid(:,10))'};

for ip = 1:size(param_info, 1)
    pname  = param_info{ip, 1};
    pcol   = param_info{ip, 2};
    levels = param_info{ip, 3};
    fprintf(rfid, '### %s\n', pname);
    fprintf(rfid, '  %-12s  %-10s  %-10s  %-10s  %-10s  %s\n', ...
        'Level', 'Mean_QC', 'Validity', 'Quality', 'Sensitiv.', 'n');
    for lev = levels(:)'
        mask = R_valid(:, pcol) == lev;
        n_lev = sum(mask);
        if n_lev == 0, continue; end
        m_qc = mean(R_valid(mask, 14), 'omitnan');
        m_v  = mean(R_valid(mask, 11), 'omitnan');
        m_q  = mean(R_valid(mask, 12), 'omitnan');
        m_s  = mean(R_valid(mask, 13), 'omitnan');
        if strcmp(pname, 'BP_HZ')
            lev_str = sprintf('%s', mat2str(BP_OPTS{lev}));
        elseif isinf(lev)
            lev_str = 'Inf';
        else
            lev_str = num2str(lev);
        end
        fprintf(rfid, '  %-12s  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %d\n', ...
            lev_str, m_qc, m_v, m_q, m_s, n_lev);
    end
    fprintf(rfid, '\n');
end

% ---- Recommendation ------------------------------------------------------
fprintf(rfid, '## Recommendation\n\n');
fprintf(rfid, 'Use the parameter combination listed under "Best Parameter\n');
fprintf(rfid, 'Combination" above for this dataset and task design.\n\n');

% Build a readable description of the best pipeline
if best_row.SIMULT
    regmode = 'simultaneous regression + bandpass filtering (regbp=2)';
else
    regmode = 'sequential bandpass-then-regression (regbp=1)';
end
if best_row.MOT24
    motmode = '24-parameter motion model (6 params + derivatives + quadratics)';
else
    motmode = '12-parameter motion model (6 params + derivatives)';
end
aroma_n = best_row.N_AROMA;
if aroma_n == 0
    aromamode = 'ICA-AROMA disabled';
elseif isinf(aroma_n)
    aromamode = 'ICA-AROMA with all components';
else
    aromamode = sprintf('ICA-AROMA with %d components', aroma_n);
end

fprintf(rfid, 'Specifically:\n');
fprintf(rfid, '  - Polynomial detrending order %d\n', best_row.POLY_ORD);
fprintf(rfid, '  - Bandpass filter %s Hz\n', best_row.BP_HZ{1});
fprintf(rfid, '  - %s\n', regmode);
fprintf(rfid, '  - %s\n', motmode);
fprintf(rfid, '  - %d aCompCor noise components\n', best_row.N_ACOMP);
fprintf(rfid, '  - %s\n\n', aromamode);

fprintf(rfid, 'To apply, set these values in run_conn_preproc.m (or a copy)\n');
fprintf(rfid, 'and re-run the full pipeline including first-level analysis.\n\n');

% ---- Caveats -------------------------------------------------------------
fprintf(rfid, '## Caveats\n\n');
n_fail = sum(cellfun(@(s) startsWith(s,'FAIL'), S));
fprintf(rfid, '- %d of %d combinations failed (%.1f%%).\n', ...
    n_fail, iter, 100*n_fail/iter);
if n_fail > 0
    fprintf(rfid, '  Failed runs are marked in the spreadsheet.  Common causes:\n');
    fprintf(rfid, '  insufficient aCompCor/AROMA components in covariate TSV files,\n');
    fprintf(rfid, '  or singular design matrices with high polynomial order.\n');
end
fprintf(rfid, '- QC scores were recomputed on CONN''s sample-voxel network\n');
fprintf(rfid, '  (~500-2000 grey-matter voxels), not the full brain volume.\n');
fprintf(rfid, '  This is standard CONN behaviour for QA plots and is\n');
fprintf(rfid, '  representative of whole-brain FC, but not identical.\n');
fprintf(rfid, '- The search used a two-stage grid: Stage 1 tested %d coarse\n', ...
    sum(R(:,2)==1));
fprintf(rfid, '  combinations; Stage 2 tested %d fine combinations around\n', ...
    sum(R(:,2)==2));
fprintf(rfid, '  the Stage 1 winner.  Interactions between the discrete\n');
fprintf(rfid, '  parameters (POLY, BP, SIMULT, MOT24) and the continuous\n');
fprintf(rfid, '  parameters (N_ACOMP, N_AROMA) are only partially captured.\n');
fprintf(rfid, '- Mean_QC is the simple arithmetic mean of Validity, Quality,\n');
fprintf(rfid, '  and Sensitivity.  All three are on [0, 1].  If one score\n');
fprintf(rfid, '  matters more for your analysis, consider weighting or\n');
fprintf(rfid, '  selecting by that score alone.\n');
fprintf(rfid, '- The DataSensitivity score depends on the number of valid\n');
fprintf(rfid, '  subjects remaining after outlier exclusion; aggressive\n');
fprintf(rfid, '  denoising may raise Validity/Quality but lower Sensitivity\n');
fprintf(rfid, '  by reducing effective DOF.\n');

fclose(rfid);
log_msg(logfid, 'Report written to %s', reportFile);

% ---- Final summary to console and log ------------------------------------
fprintf('\n');
fprintf('==============================\n');
fprintf(' OPTIMIZATION COMPLETE\n');
fprintf('==============================\n');
fprintf(' Iterations : %d\n', iter);
fprintf(' Failed     : %d\n', n_fail);
fprintf(' Best Mean_QC = %.4f\n', best_row.Mean_QC);
fprintf('   POLY=%d  BP=%s  SIM=%d  M24=%d  AC=%d  AR=%s\n', ...
    best_row.POLY_ORD, best_row.BP_HZ{1}, best_row.SIMULT, best_row.MOT24, ...
    best_row.N_ACOMP, AR_str{sort_idx(1)});
fprintf(' Results    : %s\n', xlsFile);
fprintf(' Report     : %s\n', reportFile);
fprintf(' Log        : %s\n', logFile);
fprintf('==============================\n');

log_msg(logfid, '====== OPTIMIZATION FINISHED ======');
fclose(logfid);


%% ====================== LOCAL FUNCTIONS =================================

function [s1, s2, s3, mean_qc, elapsed, status] = run_one_combo( ...
        orig_preproc, poly, bp, simult, mot24, acomp, aroma, ...
        l2_covars, sens_vars)
% RUN_ONE_COMBO  Update CONN_x.Preproc and compute three QC scores.
%
%   Modifies CONN_x.Preproc in memory (restoring from orig_preproc first),
%   then calls conn_qascores which recomputes denoising on the sample-voxel
%   network on-the-fly.  No full voxel-level denoising is run.

    global CONN_x;
    s1 = NaN; s2 = NaN; s3 = NaN; mean_qc = NaN; elapsed = 0;

    try
        t0 = tic;

        % --- Restore template and apply new parameters --------------------
        CONN_x.Preproc = orig_preproc;
        CONN_x.Preproc.detrending = poly;
        CONN_x.Preproc.filter     = bp;
        CONN_x.Preproc.regbp      = 1 + simult;   % 1=sequential, 2=simultaneous
        CONN_x.Preproc.confounds  = update_confounds( ...
            orig_preproc.confounds, mot24, acomp, aroma);

        % --- Compute scores (order matters: Validity creates QC_DOF etc.) -
        s1 = conn_qascores('DataValidity',    [], []);
        s2 = conn_qascores('DataQuality',     [], [], l2_covars, {});
        s3 = conn_qascores('DataSensitivity', [], [], sens_vars, [], 'extreme');

        if isempty(s1), s1 = NaN; end
        if isempty(s2), s2 = NaN; end
        if isempty(s3), s3 = NaN; end
        mean_qc = mean([s1, s2, s3], 'omitnan');
        elapsed = toc(t0);
        status  = 'OK';

    catch err
        elapsed = toc(t0);
        status  = ['FAIL: ', strrep(err.message, newline, ' ')];
    end
end


function conf = update_confounds(conf_orig, mot24, n_acomp, n_aroma)
% UPDATE_CONFOUNDS  Create a modified confounds struct for one iteration.
%
%   Starting from the fully-initialised original confounds struct, adjusts:
%     - realignment deriv + power (MOT24)
%     - aCompCor dimensions
%     - aroma dimensions (or removes aroma entirely when n_aroma == 0)

    conf = conf_orig;

    % ---- Realignment: derivatives + optional quadratic expansion ---------
    idx = find(strcmp(conf.names, 'realignment'), 1);
    if ~isempty(idx)
        if iscell(conf.deriv),  conf.deriv{idx} = 1;  end    % always 1st derivative
        if mot24
            if iscell(conf.power), conf.power{idx} = 2; end  % quadratic -> 24 params
        else
            if iscell(conf.power), conf.power{idx} = 1; end  % linear   -> 12 params
        end
    end

    % ---- aCompCor: number of components ----------------------------------
    idx = find(strcmp(conf.names, 'aCompCor'), 1);
    if ~isempty(idx) && iscell(conf.dimensions)
        if isempty(conf.dimensions{idx})
            conf.dimensions{idx} = n_acomp;
        else
            conf.dimensions{idx}(1) = n_acomp;
        end
    end

    % ---- ICA-AROMA -------------------------------------------------------
    if n_aroma == 0
        % Remove aroma from all confound-struct fields
        idx = find(strcmp(conf.names, 'aroma'), 1);
        if ~isempty(idx)
            flds = fieldnames(conf);
            for f = 1:numel(flds)
                val = conf.(flds{f});
                if iscell(val) && numel(val) >= idx
                    val(idx) = [];
                    conf.(flds{f}) = val;
                end
            end
        end
    else
        idx = find(strcmp(conf.names, 'aroma'), 1);
        if ~isempty(idx) && iscell(conf.dimensions)
            if isempty(conf.dimensions{idx})
                conf.dimensions{idx} = n_aroma;
            else
                conf.dimensions{idx}(1) = n_aroma;
            end
        end
    end
end


function log_msg(fid, fmt, varargin)
% LOG_MSG  Write a timestamped line to the log file and to the console.
    msg = sprintf(fmt, varargin{:});
    ts  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    line = sprintf('[%s] %s\n', ts, msg);
    fprintf('%s', line);
    if fid > 0
        fprintf(fid, '%s', line);
    end
end
