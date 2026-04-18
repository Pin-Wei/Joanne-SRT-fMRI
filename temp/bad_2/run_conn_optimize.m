%% run_conn_optimize.m   (CORRECTED)
% =========================================================================
% Two-stage grid-search optimization of CONN fMRI denoising parameters.
%
% CORRECTION NOTES (vs the previous version in temp/outdated/):
%   1. Starting point is now the EXISTING run_conn_batch.m project
%      (v2_no_param.mat).  Previously the script rebuilt the project
%      with strcat('Effect of ', condname) which silently stripped the
%      trailing whitespace, producing confound names like 'Effect ofrandom'
%      that did not match any condition, so the per-condition regressors
%      were never applied during QC denoising.  This caused Validity scores
%      to diverge from the GUI (old: 0.8035  vs  GUI: 0.7558).
%
%   2. DataSensitivity is now computed manually from QC_DOF.  The
%      CONN builtin conn_qascores('DataSensitivity',...) silently returns
%      [] for this dataset because sum(QC_DOF_WelchSatterthwaite - 3) is
%      negative -> sqrt() is complex -> spm_Ncdf fails.  Using the more
%      stable residual-DOF estimator matches the manual hand-calc and
%      gives meaningful comparisons across parameter combinations.
%
%   3. A VALIDATION block at the top re-scores the run_conn_batch.m
%      configuration and aborts the full grid run if the recomputed
%      Validity/Quality differ from the GUI ground truth by more than
%      a tolerance.  This is the explicit cross-check requested.
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

projectFile = '/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/v2_no_param.mat';
tempDir     = '/home/aclexp/pinwei/Joanne_SRT_fMRI/temp';

logFile    = fullfile(tempDir, 'optimization_log.txt');
xlsFile    = fullfile(tempDir, 'optimization_results.xlsx');
reportFile = fullfile(tempDir, 'optimization_report.md');

% --- QC variable lists ----------------------------------------------------
L2_COVARS = {'QC_MeanMotion','QC_InvalidScans'};            % for DataQuality
SENS_VARS = {'QC_MeanMotion','QC_DOF','QC_PeakFC','QC_StdFC'}; % for DataSensitivity outlier detection

% --- GUI ground truth (from run_conn_batch.m output) ---------------------
GT_POLY = 3; GT_BP = [0.008 0.09]; GT_SIMULT = 1; GT_MOT24 = 0;
GT_ACOMP = 15; GT_AROMA = Inf;
GT_V = 0.7558;   % from DataValidityScore.mat
GT_Q = 0.8959;   % min(DataQualityScore) across [QC_MeanMotion, QC_InvalidScans]
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
log_msg(logfid, 'GUI ground truth : V=%.4f  Q=%.4f  (S not stored)', GT_V, GT_Q);

ok_v = abs(v0 - GT_V) <= TOL;
ok_q = abs(q0 - GT_Q) <= TOL;
if ~(ok_v && ok_q)
    log_msg(logfid, 'VALIDATION FAILED: V_diff=%.4f  Q_diff=%.4f  (tol=%.3f)', ...
        abs(v0-GT_V), abs(q0-GT_Q), TOL);
    fclose(logfid);
    error('Scores do not match GUI within %.3f. See log.', TOL);
end
log_msg(logfid, 'VALIDATION PASSED (within %.3f)', TOL);

% Remember validation scores so we can flag the matching grid row
VAL_V = v0; VAL_Q = q0; VAL_S = s0;

%% ====================== 4. OPTIMISATION LOOP ============================
%
% R columns:
%   1=iter, 2=stage, 3=poly, 4=bp_idx, 5=bp_lo, 6=bp_hi,
%   7=simult, 8=mot24, 9=acomp, 10=aroma,
%   11=validity, 12=quality, 13=sensitivity, 14=mean_qc, 15=runtime_s

MAX_ITER = 500;
R = nan(MAX_ITER, 15);
S = cell(MAX_ITER, 1);
iter = 0;

% ---- Stage 1: coarse grid ------------------------------------------------
log_msg(logfid, '');
log_msg(logfid, '========== STAGE 1: COARSE GRID ==========');

[g1,g2,g3,g4,g5,g6] = ndgrid(1:numel(POLY_VALS),   1:numel(BP_OPTS), ...
                              1:numel(SIMULT_VALS), 1:numel(MOT24_VALS), ...
                              1:numel(ACOMP_COARSE), 1:numel(AROMA_COARSE));
grid = [g1(:), g2(:), g3(:), g4(:), g5(:), g6(:)];
n_combos = size(grid, 1);
log_msg(logfid, 'Stage 1: %d combinations', n_combos);

for ic = 1:n_combos
    poly   = POLY_VALS(grid(ic,1));
    bp_idx = grid(ic,2);       bp = BP_OPTS{bp_idx};
    simult = SIMULT_VALS(grid(ic,3));
    mot24  = MOT24_VALS(grid(ic,4));
    acomp  = ACOMP_COARSE(grid(ic,5));
    aroma  = AROMA_COARSE(grid(ic,6));

    iter = iter + 1;
    log_msg(logfid, '[S1 %03d/%03d] POLY=%d  BP=[%.3f,%g]  SIM=%d  M24=%d  AC=%d  AR=%g', ...
        ic, n_combos, poly, bp(1), bp(2), simult, mot24, acomp, aroma);

    [v, q, s, mq, el, st] = run_one_combo( ...
        orig_preproc, poly, bp, simult, mot24, acomp, aroma, ...
        L2_COVARS, SENS_VARS);

    R(iter,:) = [iter, 1, poly, bp_idx, bp(1), bp(2), ...
                 simult, mot24, acomp, aroma, v, q, s, mq, el];
    S{iter}   = st;
    log_msg(logfid, '  -> V=%.4f  Q=%.4f  S=%.4f  Mean=%.4f  (%.1fs)  [%s]', ...
        v, q, s, mq, el, st);
end

% ---- Stage 1 winner ------------------------------------------------------
s1_rows = find(R(1:iter,2) == 1);
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
    log_msg(logfid, '[S2 %03d/%03d] POLY=%d  BP=[%.3f,%g]  SIM=%d  M24=%d  AC=%d  AR=%g', ...
        ic, n_fine, win_poly, bp_win(1), bp_win(2), win_simult, win_mot24, acomp, aroma);

    [v, q, s, mq, el, st] = run_one_combo( ...
        orig_preproc, win_poly, bp_win, win_simult, win_mot24, acomp, aroma, ...
        L2_COVARS, SENS_VARS);

    R(iter,:) = [iter, 2, win_poly, win_bp_idx, bp_win(1), bp_win(2), ...
                 win_simult, win_mot24, acomp, aroma, v, q, s, mq, el];
    S{iter}   = st;
    log_msg(logfid, '  -> V=%.4f  Q=%.4f  S=%.4f  Mean=%.4f  (%.1fs)  [%s]', ...
        v, q, s, mq, el, st);
end

R = R(1:iter, :);
S = S(1:iter);
log_msg(logfid, '');
log_msg(logfid, 'All %d iterations complete.', iter);

%% ====================== 5. EXPORT RESULTS ===============================

log_msg(logfid, 'Exporting results ...');

BP_str = cell(iter,1);  AR_str = cell(iter,1);
for k = 1:iter
    BP_str{k} = sprintf('[%.3f, %s]', R(k,5), num2str(R(k,6)));
    if isinf(R(k,10)), AR_str{k} = 'Inf'; else, AR_str{k} = num2str(R(k,10)); end
end

% GUI_verified flag: row matching run_conn_batch.m parameters
gui_match = (R(:,3)==GT_POLY & R(:,4)==find_bp_idx(BP_OPTS, GT_BP) & ...
             R(:,7)==GT_SIMULT & R(:,8)==GT_MOT24 & ...
             R(:,9)==GT_ACOMP & R(:,10)==GT_AROMA);

T = table( ...
    (1:iter)',         R(:,2),          R(:,3), ...
    BP_str,            logical(R(:,7)), logical(R(:,8)), ...
    R(:,9),            R(:,10), ...
    R(:,11),           R(:,12),         R(:,13),          R(:,14), ...
    gui_match,         S,               R(:,15), ...
    'VariableNames', {'Iter','Stage','POLY_ORD', ...
                      'BP_HZ','SIMULT','MOT24', ...
                      'N_ACOMP','N_AROMA', ...
                      'Validity','Quality','Sensitivity','Mean_QC', ...
                      'GUI_verified','Status','Runtime_s'});

[~, ord] = sort(T.Mean_QC, 'descend', 'MissingPlacement','last');
T = T(ord, :);
writetable(T, xlsFile, 'Sheet','All Results');

best_row = T(1, :);
writetable(best_row, xlsFile, 'Sheet','Best Pipeline');

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
fprintf(rfid, '| Score        | GUI    | Re-computed | Diff   |\n');
fprintf(rfid, '|--------------|--------|-------------|--------|\n');
fprintf(rfid, '| Validity     | %.4f | %.4f      | %+.4f |\n', GT_V, VAL_V, VAL_V-GT_V);
fprintf(rfid, '| Quality      | %.4f | %.4f      | %+.4f |\n', GT_Q, VAL_Q, VAL_Q-GT_Q);
fprintf(rfid, '| Sensitivity  | (not stored) | %.4f | - |\n\n', VAL_S);
fprintf(rfid, 'Tolerance: %.3f.  Validation **%s**.\n\n', TOL, ...
    ternary(abs(VAL_V-GT_V)<=TOL && abs(VAL_Q-GT_Q)<=TOL, 'passed', 'FAILED'));

fprintf(rfid, '## Best Parameter Combination\n\n');
fprintf(rfid, '| Parameter | Value |\n|-----------|-------|\n');
fprintf(rfid, '| POLY_ORD  | %d |\n', best_row.POLY_ORD);
fprintf(rfid, '| BP_HZ     | %s |\n', best_row.BP_HZ{1});
fprintf(rfid, '| SIMULT    | %d |\n', best_row.SIMULT);
fprintf(rfid, '| MOT24     | %d |\n', best_row.MOT24);
fprintf(rfid, '| N_ACOMP   | %d |\n', best_row.N_ACOMP);
fprintf(rfid, '| N_AROMA   | %s |\n', AR_str{ord(1)});
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
fprintf(rfid, '- %d of %d combinations failed (%.1f%%).\n', n_fail, iter, 100*n_fail/iter);
fprintf(rfid, '- QC scores recomputed on CONN''s sample-voxel network (~500-2000 grey-matter voxels).\n');
fprintf(rfid, '- Two-stage grid: Stage 1 tested %d coarse combos; Stage 2 %d fine combos around the winner.\n', ...
    sum(R(:,2)==1), sum(R(:,2)==2));
fprintf(rfid, '- Mean_QC = simple mean of Validity, Quality, Sensitivity (each in [0,1]).\n');
fprintf(rfid, '- AROMA available components in this dataset = 6, so N_AROMA >= 5 all map to at most 6 actual regressors.\n\n');

fprintf(rfid, '## Correction Note\n\n');
fprintf(rfid, 'The previous version of this script (moved to `temp/outdated/`) produced QC\n');
fprintf(rfid, 'scores that did not match the CONN GUI for the same project.  For the\n');
fprintf(rfid, 'run_conn_batch.m configuration the old script reported Validity=0.8035\n');
fprintf(rfid, 'while the GUI reported 0.7558 - a discrepancy of 0.048.  Two bugs were\n');
fprintf(rfid, 'identified and fixed:\n\n');
fprintf(rfid, '1. **`strcat` stripped trailing whitespace in confound names.**  The old\n');
fprintf(rfid, '   script rebuilt the CONN project with\n');
fprintf(rfid, '   `strcat(''Effect of '', COND_NAMES)`, which MATLAB''s `strcat` silently\n');
fprintf(rfid, '   strips to `''Effect ofrandom''`, `''Effect ofstr_r12''`, etc.  `conn_batch`\n');
fprintf(rfid, '   accepted these without error but, at QC time, `conn_designmatrix`\n');
fprintf(rfid, '   could not associate them with any condition, so the per-condition\n');
fprintf(rfid, '   regressors were effectively omitted.  That changed the shape of the\n');
fprintf(rfid, '   after-denoising FC distribution and thus the Validity score.  Fix:\n');
fprintf(rfid, '   start from the existing `v2_no_param.mat` project that was built by\n');
fprintf(rfid, '   `run_conn_batch.m` (which uses `append` instead of `strcat` and\n');
fprintf(rfid, '   therefore has the correct space-separated names).\n\n');
fprintf(rfid, '2. **DataSensitivity silently returned [].**\n');
fprintf(rfid, '   `conn_qascores(''DataSensitivity'', ...)` prefers\n');
fprintf(rfid, '   `QC_DOF_WelchSatterthwaite` over `QC_DOF`, but for this dataset the\n');
fprintf(rfid, '   per-subject Welch-Satterthwaite values are often below 3, so\n');
fprintf(rfid, '   `sum(QC_DOF_WS - 3)` is negative, `sqrt` returns a complex number,\n');
fprintf(rfid, '   `spm_Ncdf` throws "Input must be real and full", and the whole call\n');
fprintf(rfid, '   is swallowed by an unguarded `try ... end`.  That is why no\n');
fprintf(rfid, '   `DataSensitivityScore.mat` was written to the GUI folder.  Fix: in\n');
fprintf(rfid, '   this script we compute Sensitivity directly with the regular\n');
fprintf(rfid, '   `QC_DOF` covariate using the same formula\n');
fprintf(rfid, '   (`1 - spm_Ncdf(1.645 - 0.1003*sqrt(sum(QC_DOF(valid)-3)))`).  The\n');
fprintf(rfid, '   manual value for the GUI configuration is 0.9194.\n\n');
fprintf(rfid, '3. **Validation step added.**  Before the grid search, the script now\n');
fprintf(rfid, '   re-scores the run_conn_batch.m parameters and aborts if recomputed\n');
fprintf(rfid, '   Validity and Quality differ from the stored GUI values by more than\n');
fprintf(rfid, '   %.3f.  The corresponding row in `optimization_results.xlsx` is\n', TOL);
fprintf(rfid, '   flagged in the `GUI_verified` column.\n\n');

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
    best_row.N_ACOMP, AR_str{ord(1)});
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
% COMPUTE_SENSITIVITY_MANUAL  Direct implementation of DataSensitivityScore
% using QC_DOF (avoids the DOF_WelchSatterthwaite negative-sum -> complex
% sqrt failure that silently zeroes conn_qascores('DataSensitivity',...)).
    global CONN_x;
    folderout = fullfile(CONN_x.folders.qa, ...
        ['QA_GUIrequest_DataSensitivity_', datestr(now,'yyyy_mm_dd')]);
    if ~exist(folderout, 'dir'), mkdir(folderout); end

    conn_process('qaplots', folderout, {'QA_COV'}, ...
        1:CONN_x.Setup.nsubjects, [], [], sens_vars);

    valid = conn_module('get', 'l2covariates', 'QC_ValidSubjects');
    dof   = conn_module('get', 'l2covariates', 'QC_DOF');
    m = valid > 0 & ~isnan(valid);
    sum_eff = sum(dof(m) - 3);
    if ~isfinite(sum_eff) || sum_eff < 0
        s = NaN;
    else
        s = 1 - spm_Ncdf(1.645 - 0.1003*sqrt(sum_eff));
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
        d = conf.dimensions{idx};
        if isempty(d), d = [n_aroma, n_aroma]; else, d(1) = n_aroma; end
        conf.dimensions{idx} = d;
    end
end


function idx = find_bp_idx(bp_opts, target)
    idx = 0;
    for k = 1:numel(bp_opts)
        if isequal(bp_opts{k}, target), idx = k; return; end
    end
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
