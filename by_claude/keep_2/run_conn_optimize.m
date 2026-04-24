%% run_conn_optimize.m   (CORRECTED + GAP-FILL)
% =========================================================================
% Two-stage grid-search optimization of CONN fMRI denoising parameters.
%
% CORRECTION NOTES (vs the previous version in temp/bad_2/):
%   1. Starting point is the EXISTING run_conn_batch.m project
%      (v3_no_param.mat).  The v2 version used `strcat` which silently
%      stripped trailing whitespace from confound names; v3 uses `append`
%      and the names are correctly space-separated.
%
%   2. DataSensitivity is now computed manually using the same
%      QC_DOF_WelchSatterthwaite covariate the CONN GUI uses.  The
%      previous bad_2 version used the integer QC_DOF (~108 per subject),
%      producing scores near 1.0 regardless of data quality.  The GUI
%      (conn_qascores.m line 140) prefers QC_DOF_WelchSatterthwaite
%      (fractional residual DOF, typ. 10-30 per subject) which yields
%      the 0.5421 value shown in the v3 QA output.
%
%   3. A VALIDATION block at the top re-scores the run_conn_batch.m
%      configuration and aborts the full grid run if the recomputed
%      Validity/Quality/Sensitivity differ from the GUI ground truth by
%      more than a tolerance.
%
%   4. GAP-FILL: prior run collapsed AROMA to {0} because the template
%      Preproc.confounds lacked an 'aroma' entry - but aroma IS loaded
%      as a Setup.l1covariates (run_conn_batch.m line 39), so the
%      timeseries data is available on disk.  update_confounds() now
%      ADDS an aroma entry at runtime when n_aroma>0 and aroma is
%      missing from Preproc.confounds.names.  EXISTING_RESULT_FILES
%      lets this script load prior xlsx rows, skip combinations already
%      tested, and append only new rows to a separate xlsx so the two
%      files can be compared side-by-side.
%
% Speed path (unchanged): for each parameter combination we mutate
% CONN_x.Preproc in memory and call conn_qascores, which triggers
% conn_process('qaplots',...) -> conn_qaplots.  conn_qaplots reads the
% pre-extracted COV files from disk (aCompCor/AROMA/task columns) and
% rebuilds the design matrix on-the-fly from the current CONN_x.Preproc
% (see conn_qaplots.m line 606).  No full voxel-level denoising re-run.
%
% Outputs (all in tempDir)
%   optimization_results.xlsx  - every combination tested, sorted by
%                                Mean_QC desc.  GUI_verified column flags
%                                the run_conn_batch.m-matching row.
%   optimization_report.md     - report with Correction Note section
%   optimization_log.txt       - timestamped progress log
%
% Usage
%   1. Ensure v2_no_param.mat has been fully run by run_conn_batch.m
%   2. >> run_conn_optimize
% =========================================================================

clear; clc; close all;
global CONN_x;

%% ====================== 1. CONFIGURATION ================================

projectFile = '/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/v3_no_param.mat';
tempDir     = '/home/aclexp/pinwei/Joanne_SRT_fMRI/temp';

logFile    = fullfile(tempDir, 'optimization_log.txt');
xlsFile    = fullfile(tempDir, 'compare_QC_scores_optimized_gap_fill.xlsx');
reportFile = fullfile(tempDir, 'optimization_report.md');

% --- Flexible skip logic --------------------------------------------------
% List zero or more previously-generated result files.  Each file is read,
% its parameter columns are hashed, and any combination already present is
% skipped this run.  Leave empty ({}) to test every combination.
EXISTING_RESULT_FILES = { ...
    fullfile(tempDir, 'keep', 'compare_QC_scores_optimized.xlsx')};

% --- QC variable lists ----------------------------------------------------
L2_COVARS = {'QC_MeanMotion','QC_InvalidScans'};            % for DataQuality
SENS_VARS = {'QC_MeanMotion','QC_DOF','QC_PeakFC','QC_StdFC'}; % for DataSensitivity outlier detection

% --- GUI ground truth (from run_conn_batch.m -> v3_no_param.mat) ---------
GT_POLY = 2; GT_BP = [0.008 0.09]; GT_SIMULT = 1; GT_MOT24 = 0;
GT_ACOMP = 20; GT_AROMA = 0;
GT_V = 0.7687;   % DataValidityScore.mat
GT_Q = 0.8872;   % min(DataQualityScore) = QC_MeanMotion (other is 0.9440)
GT_S = 0.5421;   % DataSensitivityScore.mat (now that GUI writes it)
TOL  = 0.01;     % absolute tolerance for validation

% --- Search grid ----------------------------------------------------------
POLY_VALS   = [1, 2, 3];
BP_OPTS     = {[0.008 0.09]; [0.008 Inf]};
SIMULT_VALS = [1, 0];
MOT24_VALS  = [1, 0];
ACOMP_ALL   = [5, 10, 15, 20, 30, 40, 50];
AROMA_ALL   = [0, 5, 10, Inf];

% Coarse Stage-1 subsets
ACOMP_COARSE = [5, 20, 50];
AROMA_COARSE = [0, Inf];

% Clean xlsx if it already exists (writetable appends otherwise)
if exist(xlsFile, 'file'), delete(xlsFile); end

%% ====================== 2. LOAD PROJECT & TEMPLATE ======================

logfid = fopen(logFile, 'w');
log_msg(logfid, '====== CONN Denoising Optimization (CORRECTED) ======');
log_msg(logfid, 'Project: %s', projectFile);

conn('load', projectFile);
assert(~isempty(CONN_x) && isfield(CONN_x,'Preproc'), 'CONN_x not loaded');

% Snapshot the fully-populated Preproc struct (correctly-named confounds)
orig_preproc = CONN_x.Preproc;
log_msg(logfid, 'Loaded: nsubjects=%d, %d confounds', ...
    CONN_x.Setup.nsubjects, numel(orig_preproc.confounds.names));
log_msg(logfid, 'Confound names: %s', ...
    strjoin(orig_preproc.confounds.names, ' | '));

N_SUBJ = CONN_x.Setup.nsubjects;

% Aroma availability check - aroma must be in Setup.l1covariates for the
% timeseries to be loadable at runtime.  If it is, we can add aroma to
% Preproc.confounds on-the-fly (see update_confounds).  If it is NOT a
% Setup covariate either, there is no aroma data anywhere in the project
% and the AROMA grid must be collapsed to {0}.
AROMA_IN_PREPROC = any(strcmp(orig_preproc.confounds.names, 'aroma'));
AROMA_IN_SETUP   = any(strcmp(CONN_x.Setup.l1covariates.names, 'aroma'));
if AROMA_IN_PREPROC
    log_msg(logfid, 'aroma is already in Preproc.confounds - will update dimensions.');
elseif AROMA_IN_SETUP
    log_msg(logfid, 'aroma is in Setup.l1covariates but not Preproc.confounds - will add at runtime.');
else
    log_msg(logfid, 'aroma absent from Setup AND Preproc -> AROMA search set to {0}.');
    AROMA_ALL    = 0;
    AROMA_COARSE = 0;
end

%% ====================== 3. VALIDATION STEP ==============================
% Re-score the exact run_conn_batch.m configuration and compare to GUI
% ground truth.  Abort if mismatch > TOL.

log_msg(logfid, '');
log_msg(logfid, '========== VALIDATION (vs GUI ground truth) ==========');

[v0, q0, s0, mq0, el0, st0] = run_one_combo( ...
    orig_preproc, GT_POLY, GT_BP, GT_SIMULT, GT_MOT24, GT_ACOMP, GT_AROMA, ...
    L2_COVARS, SENS_VARS);
log_msg(logfid, 'Recomputed scores: V=%.4f  Q=%.4f  S=%.4f  [%.1fs]', ...
    v0, q0, s0, el0);
log_msg(logfid, 'GUI ground truth : V=%.4f  Q=%.4f  S=%.4f', GT_V, GT_Q, GT_S);

ok_v = abs(v0 - GT_V) <= TOL;
ok_q = abs(q0 - GT_Q) <= TOL;
ok_s = abs(s0 - GT_S) <= TOL;
if ~(ok_v && ok_q && ok_s)
    log_msg(logfid, 'VALIDATION FAILED: V_diff=%.4f  Q_diff=%.4f  S_diff=%.4f  (tol=%.3f)', ...
        abs(v0-GT_V), abs(q0-GT_Q), abs(s0-GT_S), TOL);
    fclose(logfid);
    error('Scores do not match GUI within %.3f. See log.', TOL);
end
log_msg(logfid, 'VALIDATION PASSED (within %.3f)', TOL);

% Remember validation scores so we can flag the matching grid row
VAL_V = v0; VAL_Q = q0; VAL_S = s0;

%% ====================== 3b. LOAD PRIOR RESULTS (skip logic) =============
% Read EXISTING_RESULT_FILES, hash their parameter columns into a set.
% Also build a virtual "stage-0" block of R so Stage 1/2 winner detection
% sees the full history.

log_msg(logfid, '');
log_msg(logfid, '========== PRIOR RESULTS (skip-already-tested) ==========');

[tested_set, prev_T] = load_prior_results(EXISTING_RESULT_FILES, BP_OPTS, logfid);
N_PREV = height(prev_T);
log_msg(logfid, 'Total prior rows: %d (skip-set size: %d)', N_PREV, tested_set.Count);

%% ====================== 4. OPTIMISATION LOOP ============================
%
% R columns:
%   1=iter, 2=stage, 3=poly, 4=bp_idx, 5=bp_lo, 6=bp_hi,
%   7=simult, 8=mot24, 9=acomp, 10=aroma,
%   11=validity, 12=quality, 13=sensitivity, 14=mean_qc, 15=runtime_s

MAX_ITER = 1000;
R = nan(MAX_ITER, 15);
S = cell(MAX_ITER, 1);
iter = 0;

% Inject prior xlsx rows as virtual stage-0 entries so winner detection
% sees them.  These are filtered out at export time.
for r = 1:N_PREV
    bp_idx_p = find_bp_idx_from_str(BP_OPTS, prev_T.BP_HZ{r});
    if bp_idx_p == 0, continue; end
    bp_p     = BP_OPTS{bp_idx_p};
    aroma_p  = norm_aroma_val(prev_T.N_AROMA(r));
    iter = iter + 1;
    R(iter,:) = [iter, 0, prev_T.POLY_ORD(r), bp_idx_p, bp_p(1), bp_p(2), ...
                 double(logical(prev_T.SIMULT(r))), double(logical(prev_T.MOT24(r))), ...
                 prev_T.N_ACOMP(r), aroma_p, ...
                 prev_T.Validity(r), prev_T.Quality(r), ...
                 prev_T.Sensitivity(r), prev_T.Mean_QC(r), 0];
    S{iter} = 'PRIOR';
end
N_VIRTUAL = iter;
log_msg(logfid, 'Injected %d virtual prior-result rows into R for winner detection.', N_VIRTUAL);

% ---- Stage 1: coarse grid ------------------------------------------------
log_msg(logfid, '');
log_msg(logfid, '========== STAGE 1: COARSE GRID ==========');

[g1,g2,g3,g4,g5,g6] = ndgrid(1:numel(POLY_VALS),   1:numel(BP_OPTS), ...
                              1:numel(SIMULT_VALS), 1:numel(MOT24_VALS), ...
                              1:numel(ACOMP_COARSE), 1:numel(AROMA_COARSE));
grid = [g1(:), g2(:), g3(:), g4(:), g5(:), g6(:)];
n_combos = size(grid, 1);
log_msg(logfid, 'Stage 1: %d combinations', n_combos);

n_skipped_s1 = 0;
for ic = 1:n_combos
    poly   = POLY_VALS(grid(ic,1));
    bp_idx = grid(ic,2);       bp = BP_OPTS{bp_idx};
    simult = SIMULT_VALS(grid(ic,3));
    mot24  = MOT24_VALS(grid(ic,4));
    acomp  = ACOMP_COARSE(grid(ic,5));
    aroma  = AROMA_COARSE(grid(ic,6));

    key = make_param_key(poly, bp_idx, simult, mot24, acomp, aroma);
    if isKey(tested_set, key)
        n_skipped_s1 = n_skipped_s1 + 1;
        log_msg(logfid, '[S1 %03d/%03d] SKIP (in EXISTING)  POLY=%d  BP=[%.3f,%g]  SIM=%d  M24=%d  AC=%d  AR=%g', ...
            ic, n_combos, poly, bp(1), bp(2), simult, mot24, acomp, aroma);
        continue;
    end

    iter = iter + 1;
    log_msg(logfid, '[S1 %03d/%03d] POLY=%d  BP=[%.3f,%g]  SIM=%d  M24=%d  AC=%d  AR=%g', ...
        ic, n_combos, poly, bp(1), bp(2), simult, mot24, acomp, aroma);

    [v, q, s, mq, el, st] = run_one_combo( ...
        orig_preproc, poly, bp, simult, mot24, acomp, aroma, ...
        L2_COVARS, SENS_VARS);

    R(iter,:) = [iter, 1, poly, bp_idx, bp(1), bp(2), ...
                 simult, mot24, acomp, aroma, v, q, s, mq, el];
    S{iter}   = st;
    tested_set(key) = true;
    log_msg(logfid, '  -> V=%.4f  Q=%.4f  S=%.4f  Mean=%.4f  (%.1fs)  [%s]', ...
        v, q, s, mq, el, st);
end
log_msg(logfid, 'Stage 1 complete: %d new combos run, %d skipped as already in EXISTING.', ...
    sum(R(1:iter,2)==1), n_skipped_s1);

% ---- Stage 1 winner (across BOTH prior + new rows) ----------------------
s1_rows = find(R(1:iter,2) == 1 | R(1:iter,2) == 0);
[~, idx_win] = max(R(s1_rows, 14));
best_s1 = s1_rows(idx_win);
win_poly   = R(best_s1, 3);
win_bp_idx = R(best_s1, 4);
win_simult = R(best_s1, 7);
win_mot24  = R(best_s1, 8);
win_acomp  = R(best_s1, 9);
win_aroma  = R(best_s1, 10);

log_msg(logfid, '');
log_msg(logfid, 'Stage 1 WINNER: iter=%d  POLY=%d  BP_idx=%d  SIM=%d  M24=%d  AC=%d  AR=%g  Mean_QC=%.4f', ...
    best_s1, win_poly, win_bp_idx, win_simult, win_mot24, win_acomp, win_aroma, R(best_s1,14));

% ---- Stage 2: fine grid around winner -----------------------------------
log_msg(logfid, '');
log_msg(logfid, '========== STAGE 2: FINE GRID ==========');

acomp_fine = ACOMP_ALL(abs(ACOMP_ALL - win_acomp) <= 15);
win_pos = find(ACOMP_ALL == win_acomp, 1);
lo = max(1, win_pos - 1);  hi = min(numel(ACOMP_ALL), win_pos + 1);
acomp_fine = unique([acomp_fine, ACOMP_ALL(lo:hi)]);
aroma_fine = AROMA_ALL;

[ga, gr] = ndgrid(1:numel(acomp_fine), 1:numel(aroma_fine));
grid2 = [ga(:), gr(:)];

% Skip combos already in tested_set (covers both EXISTING xlsx and newly
% completed Stage 1 rows).
n_fine_all = size(grid2, 1);
n_skipped_s2 = 0;
log_msg(logfid, 'Stage 2: %d candidate combinations (acomp_fine=%s, aroma_fine=%s)', ...
    n_fine_all, mat2str(acomp_fine), mat2str(aroma_fine));

bp_win = BP_OPTS{win_bp_idx};

for ic = 1:n_fine_all
    acomp = acomp_fine(grid2(ic,1));
    aroma = aroma_fine(grid2(ic,2));

    key = make_param_key(win_poly, win_bp_idx, win_simult, win_mot24, acomp, aroma);
    if isKey(tested_set, key)
        n_skipped_s2 = n_skipped_s2 + 1;
        log_msg(logfid, '[S2 %03d/%03d] SKIP (already tested) AC=%d  AR=%g', ...
            ic, n_fine_all, acomp, aroma);
        continue;
    end

    iter = iter + 1;
    log_msg(logfid, '[S2 %03d/%03d] POLY=%d  BP=[%.3f,%g]  SIM=%d  M24=%d  AC=%d  AR=%g', ...
        ic, n_fine_all, win_poly, bp_win(1), bp_win(2), win_simult, win_mot24, acomp, aroma);

    [v, q, s, mq, el, st] = run_one_combo( ...
        orig_preproc, win_poly, bp_win, win_simult, win_mot24, acomp, aroma, ...
        L2_COVARS, SENS_VARS);

    R(iter,:) = [iter, 2, win_poly, win_bp_idx, bp_win(1), bp_win(2), ...
                 win_simult, win_mot24, acomp, aroma, v, q, s, mq, el];
    S{iter}   = st;
    tested_set(key) = true;
    log_msg(logfid, '  -> V=%.4f  Q=%.4f  S=%.4f  Mean=%.4f  (%.1fs)  [%s]', ...
        v, q, s, mq, el, st);
end
log_msg(logfid, 'Stage 2 complete: %d new combos run, %d skipped.', ...
    sum(R(1:iter,2)==2), n_skipped_s2);

R = R(1:iter, :);
S = S(1:iter);
n_new = sum(R(:,2) ~= 0);
n_prior = sum(R(:,2) == 0);
log_msg(logfid, '');
log_msg(logfid, 'All iterations complete: %d new + %d prior = %d total in merged set.', ...
    n_new, n_prior, iter);

%% ====================== 5. EXPORT RESULTS ===============================

log_msg(logfid, 'Exporting results ...');

BP_str_all = cell(iter,1);  AR_str_all = cell(iter,1);
for k = 1:iter
    BP_str_all{k} = sprintf('[%.3f, %s]', R(k,5), num2str(R(k,6)));
    if isinf(R(k,10)), AR_str_all{k} = 'Inf'; else, AR_str_all{k} = num2str(R(k,10)); end
end

% GUI_verified flag: row matching run_conn_batch.m parameters
gui_match_all = (R(:,3)==GT_POLY & R(:,4)==find_bp_idx(BP_OPTS, GT_BP) & ...
                 R(:,7)==GT_SIMULT & R(:,8)==GT_MOT24 & ...
                 R(:,9)==GT_ACOMP & R(:,10)==GT_AROMA);

T_all = table( ...
    (1:iter)',         R(:,2),              R(:,3), ...
    BP_str_all,        logical(R(:,7)),     logical(R(:,8)), ...
    R(:,9),            AR_str_all, ...
    R(:,11),           R(:,12),             R(:,13),          R(:,14), ...
    gui_match_all,     S,                   R(:,15), ...
    'VariableNames', {'Iter','Stage','POLY_ORD', ...
                      'BP_HZ','SIMULT','MOT24', ...
                      'N_ACOMP','N_AROMA', ...
                      'Validity','Quality','Sensitivity','Mean_QC', ...
                      'GUI_verified','Status','Runtime_s'});
% N_AROMA is stored as a cellstr ('0','5','10','Inf') because Excel's
% numeric format cannot represent Inf - MATLAB's writetable would write
% Inf as the sentinel 65535 otherwise.  Readers should str2double the
% column, mapping 'Inf' -> Inf explicitly.

% New rows only - same column structure as the prior xlsx for side-by-side
% comparison.  Stage 0 (prior) rows are excluded.
is_new_row = T_all.Stage ~= 0;
T_new = T_all(is_new_row, :);
[~, ord_new] = sort(T_new.Mean_QC, 'descend', 'MissingPlacement','last');
T_new = T_new(ord_new, :);
writetable(T_new, xlsFile, 'Sheet','All Results');

% Best Pipeline - global best across the MERGED set (new + prior)
[~, ord_all] = sort(T_all.Mean_QC, 'descend', 'MissingPlacement','last');
T_all_sorted = T_all(ord_all, :);
best_row = T_all_sorted(1, :);
best_row.Origin = string(ternary(best_row.Stage==0, 'prior', 'new'));
writetable(best_row, xlsFile, 'Sheet','Best Pipeline');

% Merged history - full picture for downstream analysis
writetable(T_all_sorted, xlsFile, 'Sheet','Merged History');

log_msg(logfid, 'Spreadsheet: %s', xlsFile);

% ---- Markdown report ----------------------------------------------------
rfid = fopen(reportFile, 'w');
fprintf(rfid, '# Denoising Pipeline Optimization Report\n\n');
fprintf(rfid, '- Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(rfid, '- Project:   `%s`\n', projectFile);
fprintf(rfid, '- Subjects:  %d\n', N_SUBJ);
fprintf(rfid, '- Iterations: %d (Stage 1: %d, Stage 2: %d)\n\n', ...
    iter, sum(R(:,2)==1), sum(R(:,2)==2));

fprintf(rfid, '## Validation Against GUI Ground Truth\n\n');
fprintf(rfid, '| Score        | GUI    | Re-computed | Diff    |\n');
fprintf(rfid, '|--------------|--------|-------------|---------|\n');
fprintf(rfid, '| Validity     | %.4f | %.4f      | %+.4f |\n', GT_V, VAL_V, VAL_V-GT_V);
fprintf(rfid, '| Quality      | %.4f | %.4f      | %+.4f |\n', GT_Q, VAL_Q, VAL_Q-GT_Q);
fprintf(rfid, '| Sensitivity  | %.4f | %.4f      | %+.4f |\n\n', GT_S, VAL_S, VAL_S-GT_S);
fprintf(rfid, 'Tolerance: %.3f.  Validation **%s**.\n\n', TOL, ...
    ternary(abs(VAL_V-GT_V)<=TOL && abs(VAL_Q-GT_Q)<=TOL && abs(VAL_S-GT_S)<=TOL, ...
            'passed', 'FAILED'));

fprintf(rfid, '## Best Parameter Combination\n\n');
fprintf(rfid, '| Parameter | Value |\n|-----------|-------|\n');
fprintf(rfid, '| POLY_ORD  | %d |\n', best_row.POLY_ORD);
fprintf(rfid, '| BP_HZ     | %s |\n', best_row.BP_HZ{1});
fprintf(rfid, '| SIMULT    | %d |\n', best_row.SIMULT);
fprintf(rfid, '| MOT24     | %d |\n', best_row.MOT24);
fprintf(rfid, '| N_ACOMP   | %d |\n', best_row.N_ACOMP);
fprintf(rfid, '| N_AROMA   | %s |\n', AR_str_all{ord_all(1)});
fprintf(rfid, '| Validity  | %.4f |\n', best_row.Validity);
fprintf(rfid, '| Quality   | %.4f |\n', best_row.Quality);
fprintf(rfid, '| Sensitivity | %.4f |\n', best_row.Sensitivity);
fprintf(rfid, '| **Mean_QC** | **%.4f** |\n\n', best_row.Mean_QC);

% Main-effect analysis
fprintf(rfid, '## Score Distributions (main-effect summary)\n\n');
fprintf(rfid, 'Mean Mean_QC across combinations grouped by each parameter level.\n\n');
valid_mask = ~isnan(R(:,14));  Rv = R(valid_mask, :);
param_info = { ...
    'POLY_ORD', 3,  POLY_VALS;
    'BP_HZ',    4,  1:numel(BP_OPTS);
    'SIMULT',   7,  SIMULT_VALS;
    'MOT24',    8,  MOT24_VALS;
    'N_ACOMP',  9,  unique(Rv(:,9))';
    'N_AROMA', 10,  unique(Rv(:,10))'};
for ip = 1:size(param_info,1)
    pname = param_info{ip,1};  pcol = param_info{ip,2};  levels = param_info{ip,3};
    fprintf(rfid, '### %s\n\n', pname);
    fprintf(rfid, '| Level | Mean_QC | Validity | Quality | Sensitivity | n |\n');
    fprintf(rfid, '|-------|---------|----------|---------|-------------|---|\n');
    for lev = levels(:)'
        mask = Rv(:,pcol) == lev;
        if sum(mask)==0, continue; end
        if strcmp(pname,'BP_HZ'),   lev_str = mat2str(BP_OPTS{lev});
        elseif isinf(lev),          lev_str = 'Inf';
        else,                       lev_str = num2str(lev);
        end
        fprintf(rfid, '| %s | %.4f | %.4f | %.4f | %.4f | %d |\n', lev_str, ...
            mean(Rv(mask,14),'omitnan'), mean(Rv(mask,11),'omitnan'), ...
            mean(Rv(mask,12),'omitnan'), mean(Rv(mask,13),'omitnan'), sum(mask));
    end
    fprintf(rfid, '\n');
end

fprintf(rfid, '## Recommendation\n\n');
fprintf(rfid, 'Use the **Best Parameter Combination** above for this dataset.\n\n');

fprintf(rfid, '## Caveats\n\n');
n_fail = sum(cellfun(@(s)startsWith(s,'FAIL'), S));
n_attempted = max(sum(R(:,2) ~= 0), 1);
fprintf(rfid, '- %d of %d new combinations failed (%.1f%%).\n', n_fail, n_attempted, 100*n_fail/n_attempted);
fprintf(rfid, '- QC scores recomputed on CONN''s sample-voxel network (~500-2000 grey-matter voxels).\n');
fprintf(rfid, '- Two-stage grid: Stage 1 tested %d coarse combos; Stage 2 %d fine combos around the winner.\n', ...
    sum(R(:,2)==1), sum(R(:,2)==2));
fprintf(rfid, '- Mean_QC = simple mean of Validity, Quality, Sensitivity (each in [0,1]).\n');
fprintf(rfid, '- AROMA available components in this dataset = 6, so N_AROMA >= 5 all map to at most 6 actual regressors.\n\n');

fprintf(rfid, '## Correction Note\n\n');
fprintf(rfid, 'This script has been corrected twice against CONN GUI ground truth.\n\n');
fprintf(rfid, '### Fix 1 (initial) - confound-name whitespace\n\n');
fprintf(rfid, 'The original `temp/outdated/` version rebuilt the CONN project via\n');
fprintf(rfid, '`strcat(''Effect of '', COND_NAMES)`.  MATLAB''s `strcat` silently strips\n');
fprintf(rfid, 'trailing whitespace, producing `''Effect ofrandom''` etc. that\n');
fprintf(rfid, '`conn_designmatrix` could not match to any condition, so the\n');
fprintf(rfid, 'per-condition regressors were never applied during QC denoising.\n');
fprintf(rfid, 'Fix: load the existing project built by `run_conn_batch.m`\n');
fprintf(rfid, '(which uses `append`, preserving the space) as the template.\n\n');
fprintf(rfid, '### Fix 2 (this pass) - DataSensitivity used the wrong DOF covariate\n\n');
fprintf(rfid, 'The `temp/bad_2/` version computed Sensitivity from `QC_DOF` (the\n');
fprintf(rfid, 'integer pre-denoising DOF, ~108 per subject).  That made\n');
fprintf(rfid, '`sum(QC_DOF(valid)-3)` ~2000 for every combination, driving the\n');
fprintf(rfid, 'score to ~0.998 regardless of the denoising parameters, and\n');
fprintf(rfid, 'disagreeing with the GUI value of **%.4f** stored in\n', GT_S);
fprintf(rfid, '`v3_no_param/results/qa/.../DataSensitivityScore.mat`.\n\n');
fprintf(rfid, 'The CONN GUI source (`conn_qascores.m` line 140-143) prefers\n');
fprintf(rfid, '`QC_DOF_WelchSatterthwaite` - the fractional residual DOF populated\n');
fprintf(rfid, 'by the `QA_COV` step (typically 10-30 per subject).  Summing\n');
fprintf(rfid, '`(QC_DOF_WS(valid) - 3)` for this dataset gives ~304.7, and\n');
fprintf(rfid, '`1 - spm_Ncdf(1.645 - 0.1003*sqrt(304.7)) = 0.5421` - matches GUI.\n\n');
fprintf(rfid, '`compute_sensitivity_manual` now pulls `QC_DOF_WelchSatterthwaite`\n');
fprintf(rfid, 'first and only falls back to `QC_DOF` if that covariate is missing.\n\n');
fprintf(rfid, '### Validation step\n\n');
fprintf(rfid, 'Before the grid search the script re-scores the run_conn_batch.m\n');
fprintf(rfid, 'configuration and aborts if recomputed Validity/Quality/Sensitivity\n');
fprintf(rfid, 'differ from the stored GUI values by more than %.3f.  The matching\n', TOL);
fprintf(rfid, 'row in `compare_QC_scores_optimized.xlsx` is flagged `GUI_verified`.\n\n');

%% --- Gap Fill section (always written; trivially short when N_PREV==0) ---
fprintf(rfid, '## Second Optimization Run - Gap Fill\n\n');

% Audit: count missing combinations from the FULL intended grid.
n_intended_s1 = numel(POLY_VALS) * numel(BP_OPTS) * numel(SIMULT_VALS) * ...
                numel(MOT24_VALS) * numel(ACOMP_COARSE) * numel(AROMA_COARSE);
n_prior_s1 = sum(R(:,2) == 0);
fprintf(rfid, '### Audit Results\n\n');
fprintf(rfid, '- Prior xlsx files merged: %d (rows: %d).\n', ...
    numel(EXISTING_RESULT_FILES), n_prior_s1);
fprintf(rfid, '- Stage 1 intended (full coarse grid): %d combinations.\n', n_intended_s1);
fprintf(rfid, '- Stage 1 newly run this pass: %d combinations.\n', sum(R(:,2)==1));
fprintf(rfid, '- Stage 2 newly run this pass: %d combinations.\n', sum(R(:,2)==2));
fprintf(rfid, '- Total new combinations: %d.\n\n', n_attempted);

fprintf(rfid, '### Root Cause of Missing N_AROMA > 0 Runs\n\n');
fprintf(rfid, 'The previous run collapsed AROMA -> {0} because\n');
fprintf(rfid, '`Preproc.confounds.names` lacked an `aroma` entry (the v3 project was\n');
fprintf(rfid, 'built by `run_conn_batch.m` with `N_AROMA = 0`, so line 47 of that\n');
fprintf(rfid, 'script never added `aroma` to `CONFOUND_NAMES`).  The auto-collapse\n');
fprintf(rfid, 'guard at line 100-106 of the prior `run_conn_optimize.m` saw the\n');
fprintf(rfid, 'absence and disabled all AROMA>0 testing.\n\n');
fprintf(rfid, 'However `aroma` IS registered in `CONN_x.Setup.l1covariates.names`\n');
fprintf(rfid, '(index 5), with all 20 subjects x 8 sessions of TSV timeseries\n');
fprintf(rfid, 'already loaded by the run_conn_batch.m Setup step.  The fix:\n');
fprintf(rfid, '`update_confounds()` now ADDS an `aroma` entry to\n');
fprintf(rfid, '`Preproc.confounds` at runtime when `n_aroma>0` and aroma is missing\n');
fprintf(rfid, 'from the names list.  The aroma availability check now consults\n');
fprintf(rfid, '`Setup.l1covariates.names` instead of `Preproc.confounds.names`.\n\n');

fprintf(rfid, '### Skip Logic\n\n');
fprintf(rfid, 'A new `EXISTING_RESULT_FILES` cell array at the top of the script\n');
fprintf(rfid, 'lists prior xlsx files.  At startup their parameter columns are\n');
fprintf(rfid, 'hashed into a set; combinations already in the set are logged and\n');
fprintf(rfid, 'skipped.  This run loaded %d prior rows from %d file(s) and skipped\n', ...
    n_prior_s1, numel(EXISTING_RESULT_FILES));
fprintf(rfid, 'the corresponding combinations during the coarse grid pass.\n\n');

fprintf(rfid, '### New Best Parameter Combination (merged set)\n\n');
fprintf(rfid, '| Parameter | Value |\n|-----------|-------|\n');
fprintf(rfid, '| Origin    | %s |\n', char(best_row.Origin));
fprintf(rfid, '| POLY_ORD  | %d |\n', best_row.POLY_ORD);
fprintf(rfid, '| BP_HZ     | %s |\n', best_row.BP_HZ{1});
fprintf(rfid, '| SIMULT    | %d |\n', best_row.SIMULT);
fprintf(rfid, '| MOT24     | %d |\n', best_row.MOT24);
fprintf(rfid, '| N_ACOMP   | %d |\n', best_row.N_ACOMP);
fprintf(rfid, '| N_AROMA   | %s |\n', AR_str_all{ord_all(1)});
fprintf(rfid, '| Validity  | %.4f |\n', best_row.Validity);
fprintf(rfid, '| Quality   | %.4f |\n', best_row.Quality);
fprintf(rfid, '| Sensitivity | %.4f |\n', best_row.Sensitivity);
fprintf(rfid, '| **Mean_QC** | **%.4f** |\n\n', best_row.Mean_QC);

fprintf(rfid, '### Comparison: Top 10 from Prior vs. New Top 10\n\n');
T_prior = T_all_sorted(T_all_sorted.Stage==0, :);
T_new10 = T_new(1:min(10,height(T_new)), :);
T_pri10 = T_prior(1:min(10,height(T_prior)), :);
fprintf(rfid, '#### Prior Top 10\n\n');
fprintf(rfid, '| POLY | BP | SIMULT | MOT24 | N_ACOMP | N_AROMA | Validity | Quality | Sens | Mean_QC |\n');
fprintf(rfid, '|------|----|--------|-------|---------|---------|----------|---------|------|---------|\n');
for r = 1:height(T_pri10)
    ar_label = T_pri10.N_AROMA{r};
    fprintf(rfid, '| %d | %s | %d | %d | %d | %s | %.4f | %.4f | %.4f | %.4f |\n', ...
        T_pri10.POLY_ORD(r), T_pri10.BP_HZ{r}, T_pri10.SIMULT(r), T_pri10.MOT24(r), ...
        T_pri10.N_ACOMP(r), ar_label, T_pri10.Validity(r), T_pri10.Quality(r), ...
        T_pri10.Sensitivity(r), T_pri10.Mean_QC(r));
end
fprintf(rfid, '\n#### New Top 10\n\n');
fprintf(rfid, '| POLY | BP | SIMULT | MOT24 | N_ACOMP | N_AROMA | Validity | Quality | Sens | Mean_QC |\n');
fprintf(rfid, '|------|----|--------|-------|---------|---------|----------|---------|------|---------|\n');
for r = 1:height(T_new10)
    ar_label = T_new10.N_AROMA{r};
    fprintf(rfid, '| %d | %s | %d | %d | %d | %s | %.4f | %.4f | %.4f | %.4f |\n', ...
        T_new10.POLY_ORD(r), T_new10.BP_HZ{r}, T_new10.SIMULT(r), T_new10.MOT24(r), ...
        T_new10.N_ACOMP(r), ar_label, T_new10.Validity(r), T_new10.Quality(r), ...
        T_new10.Sensitivity(r), T_new10.Mean_QC(r));
end
fprintf(rfid, '\n');

fclose(rfid);
log_msg(logfid, 'Report: %s', reportFile);

% ---- Console summary ----------------------------------------------------
fprintf('\n==============================\n');
fprintf(' OPTIMIZATION COMPLETE\n');
fprintf('==============================\n');
fprintf(' Iterations : %d  (failed: %d)\n', iter, n_fail);
fprintf(' Best Mean_QC = %.4f\n', best_row.Mean_QC);
fprintf('   POLY=%d  BP=%s  SIM=%d  M24=%d  AC=%d  AR=%s\n', ...
    best_row.POLY_ORD, best_row.BP_HZ{1}, best_row.SIMULT, best_row.MOT24, ...
    best_row.N_ACOMP, AR_str_all{ord_all(1)});
fprintf(' XLSX   : %s\n', xlsFile);
fprintf(' Report : %s\n', reportFile);
fprintf(' Log    : %s\n', logFile);
fprintf('==============================\n');

log_msg(logfid, '====== OPTIMIZATION FINISHED ======');
fclose(logfid);


%% ====================== LOCAL FUNCTIONS =================================

function [v, q, s, mean_qc, elapsed, status] = run_one_combo( ...
        orig_preproc, poly, bp, simult, mot24, acomp, aroma, ...
        l2_covars, sens_vars)
% RUN_ONE_COMBO  Update CONN_x.Preproc and compute the three QC scores.
%
%   Mutates CONN_x.Preproc in memory (restoring from orig_preproc first),
%   then calls conn_qascores for Validity and Quality.  Sensitivity is
%   computed manually from QC_DOF to avoid the spm_Ncdf complex-sqrt
%   issue that silently zeroes the builtin.
    global CONN_x;
    v = NaN; q = NaN; s = NaN; mean_qc = NaN; elapsed = 0;

    try
        t0 = tic;

        CONN_x.Preproc = orig_preproc;
        CONN_x.Preproc.detrending = poly;
        CONN_x.Preproc.filter     = bp;
        CONN_x.Preproc.regbp      = 1 + simult;  % 1=sequential, 2=simultaneous
        CONN_x.Preproc.confounds  = update_confounds( ...
            orig_preproc.confounds, mot24, acomp, aroma);

        v = safe_score(@() conn_qascores('DataValidity', [], []));
        q = safe_score(@() conn_qascores('DataQuality',  [], [], l2_covars, {}));
        s = compute_sensitivity_manual(sens_vars);

        mean_qc = mean([v, q, s], 'omitnan');
        elapsed = toc(t0);
        status  = 'OK';

    catch err
        elapsed = toc(t0);
        status  = ['FAIL: ', strrep(err.message, newline, ' ')];
    end
end


function out = safe_score(fn)
% SAFE_SCORE  Call a conn_qascores handle and normalise [] -> NaN.
    out = fn();
    if isempty(out), out = NaN; end
end


function s = compute_sensitivity_manual(sens_vars)
% COMPUTE_SENSITIVITY_MANUAL  Direct replication of the CONN GUI
% DataSensitivityScore (conn_qascores.m case 'datasensitivity', lines
% 118-150).  The GUI calls conn_process('qaplots', ..., {'QA_COV'}, ...)
% which populates three l2covariates for us:
%   - QC_OutlierScore          (|z| across the four SENS_VARS)
%   - QC_ValidSubjects         (1 if QC_OutlierScore <= 3, else 0)
%   - QC_DOF_WelchSatterthwaite  (residual DOF after denoising; fractional)
%
% Then applies:
%   DataSensitivityScore = 1 - spm_Ncdf(1.645 - 0.1003 * sqrt(sum(DOF_WS(valid)-3)))
%
% bad_2 wrongly used QC_DOF (the integer pre-denoising DOF, ~108 per
% subject) which makes sum_eff huge and drives the score to ~1.0 for
% every combination.
    global CONN_x;
    folderout = fullfile(CONN_x.folders.qa, ...
        ['QA_GUIrequest_DataSensitivity_', datestr(now,'yyyy_mm_dd')]);
    if ~exist(folderout, 'dir'), mkdir(folderout); end

    conn_process('qaplots', folderout, {'QA_COV'}, ...
        1:CONN_x.Setup.nsubjects, [], [], sens_vars);

    valid = conn_module('get', 'l2covariates', 'QC_ValidSubjects');
    dof   = conn_module('get', 'l2covariates', 'QC_DOF_WelchSatterthwaite');
    if isempty(dof)
        dof = conn_module('get', 'l2covariates', 'QC_DOF');  % fallback only
    end
    m = valid > 0 & ~isnan(valid);
    sum_eff = sum(dof(m) - 3);
    if ~isfinite(sum_eff) || sum_eff < 0
        s = NaN;                                  % matches GUI silent-fail
    else
        s = 1 - spm_Ncdf(1.645 - 0.1003 * sqrt(sum_eff));
    end
end


function conf = update_confounds(conf_orig, mot24, n_acomp, n_aroma)
% UPDATE_CONFOUNDS  Derive a confounds struct for one iteration.
%   Assumes conf_orig was extracted from a fully-denoised project so the
%   'realignment', 'aCompCor' and 'aroma' entries are all present and
%   'dimensions{idx}' is a 1x2 vector [requested, available].
    conf = conf_orig;

    % Realignment: linear vs quadratic expansion (12 vs 24 parameters)
    idx = find(strcmp(conf.names, 'realignment'), 1);
    if ~isempty(idx)
        if iscell(conf.deriv),  conf.deriv{idx} = 1;  end
        if iscell(conf.power),  conf.power{idx} = 1 + double(logical(mot24)); end
    end

    % aCompCor components
    idx = find(strcmp(conf.names, 'aCompCor'), 1);
    if ~isempty(idx) && iscell(conf.dimensions)
        d = conf.dimensions{idx};
        if isempty(d), d = [n_acomp, n_acomp]; else, d(1) = n_acomp; end
        conf.dimensions{idx} = d;
    end

    % ICA-AROMA
    idx = find(strcmp(conf.names, 'aroma'), 1);
    if n_aroma == 0
        % Remove aroma entirely from every cell-array field
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
    elseif ~isempty(idx) && iscell(conf.dimensions)
        % aroma already present - just update the requested dimension
        d = conf.dimensions{idx};
        if isempty(d), d = [n_aroma, n_aroma]; else, d(1) = n_aroma; end
        conf.dimensions{idx} = d;
    elseif isempty(idx)
        % aroma absent from Preproc.confounds but caller wants n_aroma>0.
        % Append a new entry; CONN will pull the timeseries from
        % Setup.l1covariates at design-matrix time.
        conf.names{end+1}      = 'aroma';
        conf.types{end+1}      = 'cov';
        conf.power{end+1}      = 1;
        conf.deriv{end+1}      = false;
        conf.dimensions{end+1} = [n_aroma, n_aroma];
        conf.filter{end+1}     = 0;
        conf.fixed{end+1}      = 0;
    end
end


function idx = find_bp_idx(bp_opts, target)
    idx = 0;
    for k = 1:numel(bp_opts)
        if isequal(bp_opts{k}, target), idx = k; return; end
    end
end


function idx = find_bp_idx_from_str(bp_opts, s)
% FIND_BP_IDX_FROM_STR  Match a stored BP_HZ string back to a BP_OPTS index.
    idx = 0;
    if iscell(s), s = s{1}; end
    s = strtrim(char(s));
    for k = 1:numel(bp_opts)
        ref = sprintf('[%.3f, %s]', bp_opts{k}(1), num2str(bp_opts{k}(2)));
        if strcmp(s, ref), idx = k; return; end
    end
end


function v = norm_aroma_val(x)
% NORM_AROMA_VAL  Coerce an N_AROMA cell into the canonical numeric form.
    if iscell(x), x = x{1}; end
    if isnumeric(x)
        v = double(x);
    elseif ischar(x) || isstring(x)
        s = strtrim(char(x));
        if strcmpi(s, 'inf'), v = Inf;
        else, v = str2double(s); end
    else
        v = NaN;
    end
end


function k = make_param_key(poly, bp_idx, simult, mot24, acomp, aroma)
% MAKE_PARAM_KEY  Canonical hashable string identifying a parameter combo.
    if isinf(aroma), aroma_str = 'Inf'; else, aroma_str = sprintf('%d', round(aroma)); end
    k = sprintf('%d|%d|%d|%d|%d|%s', round(poly), round(bp_idx), ...
                double(logical(simult)), double(logical(mot24)), ...
                round(acomp), aroma_str);
end


function [tested_set, prev_T] = load_prior_results(file_list, BP_OPTS, logfid)
% LOAD_PRIOR_RESULTS  Read EXISTING_RESULT_FILES, build a tested-set hash
% and a vertcat'd table of all prior rows for virtual-row injection.
    tested_set = containers.Map('KeyType','char','ValueType','logical');
    prev_T = table();
    if isempty(file_list)
        log_msg(logfid, 'EXISTING_RESULT_FILES empty - no skip set built.');
        return;
    end
    for iFile = 1:numel(file_list)
        fpath = file_list{iFile};
        if ~exist(fpath, 'file')
            log_msg(logfid, '  WARN  file not found: %s', fpath);
            continue;
        end
        try
            Tprev = readtable(fpath, 'Sheet','All Results');
        catch
            try
                Tprev = readtable(fpath);
            catch err
                log_msg(logfid, '  WARN  failed to read %s: %s', fpath, err.message);
                continue;
            end
        end
        log_msg(logfid, '  loaded %d rows from %s', height(Tprev), fpath);
        for r = 1:height(Tprev)
            bp_idx = find_bp_idx_from_str(BP_OPTS, Tprev.BP_HZ{r});
            aroma_v = norm_aroma_val(Tprev.N_AROMA(r));
            key = make_param_key(Tprev.POLY_ORD(r), bp_idx, ...
                                 Tprev.SIMULT(r), Tprev.MOT24(r), ...
                                 Tprev.N_ACOMP(r), aroma_v);
            tested_set(key) = true;
        end
        prev_T = vertcat_safe(prev_T, Tprev);
    end
end


function T = vertcat_safe(A, B)
% VERTCAT_SAFE  Concatenate two tables tolerating column-set differences
% by aligning to the intersection of variable names.
    if isempty(A), T = B; return; end
    if isempty(B), T = A; return; end
    common = intersect(A.Properties.VariableNames, B.Properties.VariableNames, 'stable');
    T = [A(:, common); B(:, common)];
end


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end


function log_msg(fid, fmt, varargin)
    msg  = sprintf(fmt, varargin{:});
    ts   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    line = sprintf('[%s] %s\n', ts, msg);
    fprintf('%s', line);
    if fid > 0, fprintf(fid, '%s', line); end
end
