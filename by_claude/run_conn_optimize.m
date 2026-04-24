%% run_conn_optimize.m    (adds RM_GMR / GSR to the search grid)
% =========================================================================
% Rewrites by_claude/keep_2/run_conn_optimize.m for the current project
% (no_PM_260423.mat) and adds a 7th parameter:
%
%   RM_GMR  -- binary flag, 1 = add 'Grey Matter' as a Preproc.confounds
%              entry (dim=1, deriv=0, power=1). This is CONN's equivalent
%              of Global Signal Regression (GSR), motivated by Ciric 2017 /
%              Mascali 2021 benchmarks showing GSR equalises task/rest
%              motion artifacts and strengthens FC-behaviour associations.
%
% Reference project     : conn_out/no_PM_260423.mat (current, Apr 23 run;
%                         setup + denoising + QA complete; first-level
%                         aborted due to unrelated COND_LIST bug -- not
%                         needed for QC optimisation).
% Reference old script  : by_claude/keep_2/run_conn_optimize.m
% Reference pipeline cfg: scripts/run_conn_batch.m
%
% Design preserved from keep_2 version:
%   - Mutates CONN_x.Preproc in memory and calls conn_qascores (fast path;
%     no voxel-level re-denoising).
%   - DataSensitivity via QC_DOF_WelchSatterthwaite (GUI-matching path).
%   - Validates against current project's stored GUI scores before the
%     grid search; aborts on mismatch > TOL.
%   - Exports xlsx (All Results / Best Pipeline / Merged History sheets)
%     plus markdown report and timestamped log.
%
% Note on skip-set: the keep_2 xlsx was produced against v3_no_param.mat,
% a pipeline WITHOUT aroma and without GSR. Its Mean_QC values are not
% directly comparable to the current project, so EXISTING_RESULT_FILES is
% intentionally empty -- all 288 coarse combinations are run fresh.
%
% Usage
%   matlab -batch "run('/home/aclexp/pinwei/Joanne_SRT_fMRI/scripts/run_conn_optimize.m')"
% =========================================================================

clear; clc; close all;
global CONN_x;

%% ====================== 1. CONFIGURATION ================================

projectFile = '/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/no_PM_260423.mat';
tempDir     = '/home/aclexp/pinwei/Joanne_SRT_fMRI/temp';
if ~isfolder(tempDir), mkdir(tempDir); end

logFile    = fullfile(tempDir, 'optimization_log.txt');
xlsFile    = fullfile(tempDir, 'optimization_results.xlsx');
reportFile = fullfile(tempDir, 'optimization_report.md');

EXISTING_RESULT_FILES = {};   % see note in header

L2_COVARS = {'QC_MeanMotion','QC_InvalidScans'};
SENS_VARS = {'QC_MeanMotion','QC_DOF','QC_PeakFC','QC_StdFC'};

% --- GUI ground truth (no_PM_260423.mat baseline) -----------------------
% Matches the v5 pipeline: run_conn_batch.m with MOT24=1, SIMULT=0,
% POLY=3, BP=[0.008 0.09], ACOMP=20, AROMA=Inf, no GSR.
GT_POLY = 3; GT_BP = [0.008 0.09]; GT_SIMULT = 0; GT_MOT24 = 1;
GT_ACOMP = 20; GT_AROMA = Inf; GT_RMGMR = 0;
GT_V = 0.7133; GT_Q = 0.9120; GT_S = 1.0000;
TOL  = 0.01;

% --- Search grid --------------------------------------------------------
% Full parameter space (for reporting); stage-1 uses the coarse subsets
% below to keep runtime tractable (~60s per combo when AROMA=Inf).
POLY_VALS    = [1, 2, 3];
BP_OPTS      = {[0.008 0.09]; [0.008 Inf]};
SIMULT_VALS  = [1, 0];
MOT24_VALS   = [1, 0];
ACOMP_ALL    = [5, 10, 15, 20, 30, 40, 50];
AROMA_ALL    = [0, 5, 10, Inf];
RM_GMR_VALS  = [0, 1];                                 % <-- NEW

% Coarse subsets. POLY_COARSE drops POLY=1 (Mean_QC=0.80 in the v3 old
% run vs POLY=3's 0.81) to shave ~1/3 of stage-1 combos; stage 2 stays
% flexible.  Result: 2*2*2*2*3*2*2 = 192 combos.
POLY_COARSE   = [2, 3];
ACOMP_COARSE  = [5, 20, 50];
AROMA_COARSE  = [0, Inf];
RM_GMR_COARSE = [0, 1];

if exist(xlsFile, 'file'), delete(xlsFile); end

%% ====================== 2. LOAD PROJECT & TEMPLATE ======================

logfid = fopen(logFile, 'w');
log_msg(logfid, '====== CONN Denoising Optimization (with RM_GMR / GSR) ======');
log_msg(logfid, 'Project: %s', projectFile);

conn('load', projectFile);
assert(~isempty(CONN_x) && isfield(CONN_x,'Preproc'), 'CONN_x not loaded');

orig_preproc = CONN_x.Preproc;
log_msg(logfid, 'Loaded: nsubjects=%d, %d confounds', ...
    CONN_x.Setup.nsubjects, numel(orig_preproc.confounds.names));
log_msg(logfid, 'Confound names: %s', strjoin(orig_preproc.confounds.names, ' | '));
N_SUBJ = CONN_x.Setup.nsubjects;

AROMA_IN_PREPROC = any(strcmp(orig_preproc.confounds.names, 'aroma'));
AROMA_IN_SETUP   = any(strcmp(CONN_x.Setup.l1covariates.names, 'aroma'));
if AROMA_IN_PREPROC
    log_msg(logfid, 'aroma is in Preproc.confounds -> will update dimensions.');
elseif AROMA_IN_SETUP
    log_msg(logfid, 'aroma is in Setup.l1covariates -> will add to Preproc at runtime.');
else
    log_msg(logfid, 'aroma absent from Setup AND Preproc -> AROMA search collapsed to {0}.');
    AROMA_ALL = 0; AROMA_COARSE = 0;
end

GM_IN_SETUP = any(strcmp(CONN_x.Setup.rois.names, 'Grey Matter'));
if ~GM_IN_SETUP
    log_msg(logfid, 'Grey Matter ROI absent from Setup.rois -> RM_GMR forced to 0.');
    RM_GMR_VALS = 0; RM_GMR_COARSE = 0;
else
    log_msg(logfid, 'Grey Matter ROI present in Setup.rois -> RM_GMR grid = [0 1].');
end

%% ====================== 3. VALIDATION ==================================

log_msg(logfid, ''); log_msg(logfid, '========== VALIDATION (vs GUI ground truth) ==========');
[v0, q0, s0, ~, el0, ~] = run_one_combo( ...
    orig_preproc, GT_POLY, GT_BP, GT_SIMULT, GT_MOT24, GT_ACOMP, GT_AROMA, GT_RMGMR, ...
    L2_COVARS, SENS_VARS);
log_msg(logfid, 'Recomputed: V=%.4f  Q=%.4f  S=%.4f  [%.1fs]', v0,q0,s0,el0);
log_msg(logfid, 'GUI truth : V=%.4f  Q=%.4f  S=%.4f', GT_V,GT_Q,GT_S);

ok_v = abs(v0 - GT_V) <= TOL;
ok_q = abs(q0 - GT_Q) <= TOL;
ok_s = abs(s0 - GT_S) <= TOL;
if ~(ok_v && ok_q && ok_s)
    log_msg(logfid, 'VALIDATION FAILED: V_d=%.4f  Q_d=%.4f  S_d=%.4f  (tol=%.3f)', ...
        abs(v0-GT_V), abs(q0-GT_Q), abs(s0-GT_S), TOL);
    fclose(logfid);
    error('Scores do not match GUI within %.3f. See log.', TOL);
end
log_msg(logfid, 'VALIDATION PASSED (within %.3f)', TOL);
VAL_V = v0; VAL_Q = q0; VAL_S = s0;

%% ====================== 3b. PRIOR RESULTS (skip) =======================

log_msg(logfid, ''); log_msg(logfid, '========== PRIOR RESULTS ==========');
[tested_set, prev_T] = load_prior_results(EXISTING_RESULT_FILES, BP_OPTS, logfid);
N_PREV = height(prev_T);
log_msg(logfid, 'Prior rows: %d (skip-set: %d)', N_PREV, tested_set.Count);

%% ====================== 4. OPTIMISATION LOOP ===========================
% R columns (16):
%  1=iter 2=stage 3=poly 4=bp_idx 5=bp_lo 6=bp_hi
%  7=simult 8=mot24 9=acomp 10=aroma 11=rm_gmr
%  12=validity 13=quality 14=sensitivity 15=mean_qc 16=runtime_s

MAX_ITER = 2000;
R = nan(MAX_ITER, 16);
S = cell(MAX_ITER, 1);
iter = 0;

for r = 1:N_PREV
    bp_idx_p = find_bp_idx_from_str(BP_OPTS, prev_T.BP_HZ{r});
    if bp_idx_p == 0, continue; end
    bp_p    = BP_OPTS{bp_idx_p};
    aroma_p = norm_aroma_val(prev_T.N_AROMA(r));
    iter = iter + 1;
    R(iter,:) = [iter, 0, prev_T.POLY_ORD(r), bp_idx_p, bp_p(1), bp_p(2), ...
                 double(logical(prev_T.SIMULT(r))), double(logical(prev_T.MOT24(r))), ...
                 prev_T.N_ACOMP(r), aroma_p, 0, ...
                 prev_T.Validity(r), prev_T.Quality(r), prev_T.Sensitivity(r), ...
                 prev_T.Mean_QC(r), 0];
    S{iter} = 'PRIOR';
end
N_VIRTUAL = iter;
log_msg(logfid, 'Injected %d prior rows as virtual stage-0 entries.', N_VIRTUAL);

% ---- Stage 1 coarse grid ----------------------------------------------
log_msg(logfid, ''); log_msg(logfid, '========== STAGE 1: COARSE GRID ==========');
[g1,g2,g3,g4,g5,g6,g7] = ndgrid(1:numel(POLY_COARSE),  1:numel(BP_OPTS), ...
                                 1:numel(SIMULT_VALS),  1:numel(MOT24_VALS), ...
                                 1:numel(ACOMP_COARSE), 1:numel(AROMA_COARSE), ...
                                 1:numel(RM_GMR_COARSE));
grid = [g1(:), g2(:), g3(:), g4(:), g5(:), g6(:), g7(:)];
n_combos = size(grid, 1);
log_msg(logfid, 'Stage 1: %d combinations', n_combos);

n_skipped_s1 = 0;
for ic = 1:n_combos
    poly   = POLY_COARSE(grid(ic,1));
    bp_idx = grid(ic,2);  bp = BP_OPTS{bp_idx};
    simult = SIMULT_VALS(grid(ic,3));
    mot24  = MOT24_VALS(grid(ic,4));
    acomp  = ACOMP_COARSE(grid(ic,5));
    aroma  = AROMA_COARSE(grid(ic,6));
    rm_gmr = RM_GMR_COARSE(grid(ic,7));

    key = make_param_key(poly, bp_idx, simult, mot24, acomp, aroma, rm_gmr);
    if isKey(tested_set, key)
        n_skipped_s1 = n_skipped_s1 + 1;
        log_msg(logfid, '[S1 %03d/%03d] SKIP  POLY=%d  BP=[%.3f,%g]  SIM=%d  M24=%d  AC=%d  AR=%g  GMR=%d', ...
            ic,n_combos,poly,bp(1),bp(2),simult,mot24,acomp,aroma,rm_gmr);
        continue;
    end
    iter = iter + 1;
    log_msg(logfid, '[S1 %03d/%03d] POLY=%d  BP=[%.3f,%g]  SIM=%d  M24=%d  AC=%d  AR=%g  GMR=%d', ...
        ic,n_combos,poly,bp(1),bp(2),simult,mot24,acomp,aroma,rm_gmr);

    [v,q,s,mq,el,st] = run_one_combo( ...
        orig_preproc, poly, bp, simult, mot24, acomp, aroma, rm_gmr, L2_COVARS, SENS_VARS);
    R(iter,:) = [iter, 1, poly, bp_idx, bp(1), bp(2), simult, mot24, acomp, aroma, rm_gmr, v, q, s, mq, el];
    S{iter}   = st;
    tested_set(key) = true;
    log_msg(logfid, '  -> V=%.4f  Q=%.4f  S=%.4f  Mean=%.4f  (%.1fs)  [%s]', v,q,s,mq,el,st);
end
log_msg(logfid, 'Stage 1: %d new, %d skipped.', sum(R(1:iter,2)==1), n_skipped_s1);

s1_rows = find(R(1:iter,2) == 1 | R(1:iter,2) == 0);
[~, idx_win] = max(R(s1_rows, 15));
best_s1 = s1_rows(idx_win);
win_poly   = R(best_s1, 3);
win_bp_idx = R(best_s1, 4);
win_simult = R(best_s1, 7);
win_mot24  = R(best_s1, 8);
win_acomp  = R(best_s1, 9);
win_aroma  = R(best_s1, 10);
win_rm_gmr = R(best_s1, 11);
log_msg(logfid, '');
log_msg(logfid, 'Stage 1 WINNER: POLY=%d  BP_idx=%d  SIM=%d  M24=%d  AC=%d  AR=%g  GMR=%d  Mean=%.4f', ...
    win_poly, win_bp_idx, win_simult, win_mot24, win_acomp, win_aroma, win_rm_gmr, R(best_s1,15));

% ---- Stage 2 fine grid ------------------------------------------------
log_msg(logfid, ''); log_msg(logfid, '========== STAGE 2: FINE GRID ==========');
acomp_fine = ACOMP_ALL(abs(ACOMP_ALL - win_acomp) <= 15);
win_pos = find(ACOMP_ALL == win_acomp, 1);
if ~isempty(win_pos)
    lo = max(1, win_pos - 1); hi = min(numel(ACOMP_ALL), win_pos + 1);
    acomp_fine = unique([acomp_fine, ACOMP_ALL(lo:hi)]);
end
aroma_fine  = AROMA_ALL;
rm_gmr_fine = RM_GMR_VALS;
[ga, gr, gm] = ndgrid(1:numel(acomp_fine), 1:numel(aroma_fine), 1:numel(rm_gmr_fine));
grid2 = [ga(:), gr(:), gm(:)];
n_fine_all = size(grid2, 1);
log_msg(logfid, 'Stage 2: %d candidate combos (AC=%s  AR=%s  GMR=%s)', ...
    n_fine_all, mat2str(acomp_fine), mat2str(aroma_fine), mat2str(rm_gmr_fine));
bp_win = BP_OPTS{win_bp_idx};
n_skipped_s2 = 0;
for ic = 1:n_fine_all
    acomp  = acomp_fine(grid2(ic,1));
    aroma  = aroma_fine(grid2(ic,2));
    rm_gmr = rm_gmr_fine(grid2(ic,3));
    key = make_param_key(win_poly, win_bp_idx, win_simult, win_mot24, acomp, aroma, rm_gmr);
    if isKey(tested_set, key)
        n_skipped_s2 = n_skipped_s2 + 1;
        log_msg(logfid, '[S2 %03d/%03d] SKIP  AC=%d  AR=%g  GMR=%d', ic,n_fine_all,acomp,aroma,rm_gmr);
        continue;
    end
    iter = iter + 1;
    log_msg(logfid, '[S2 %03d/%03d] AC=%d  AR=%g  GMR=%d', ic,n_fine_all,acomp,aroma,rm_gmr);
    [v,q,s,mq,el,st] = run_one_combo( ...
        orig_preproc, win_poly, bp_win, win_simult, win_mot24, acomp, aroma, rm_gmr, L2_COVARS, SENS_VARS);
    R(iter,:) = [iter, 2, win_poly, win_bp_idx, bp_win(1), bp_win(2), ...
                 win_simult, win_mot24, acomp, aroma, rm_gmr, v, q, s, mq, el];
    S{iter}   = st;
    tested_set(key) = true;
    log_msg(logfid, '  -> V=%.4f  Q=%.4f  S=%.4f  Mean=%.4f  (%.1fs)  [%s]', v,q,s,mq,el,st);
end
log_msg(logfid, 'Stage 2: %d new, %d skipped.', sum(R(1:iter,2)==2), n_skipped_s2);

R = R(1:iter, :); S = S(1:iter);
n_new = sum(R(:,2) ~= 0);
n_prior = sum(R(:,2) == 0);
log_msg(logfid, ''); log_msg(logfid, 'All iterations: %d new + %d prior = %d total.', n_new,n_prior,iter);

%% ====================== 5. EXPORT ======================================

log_msg(logfid, 'Exporting ...');
BP_str_all = cell(iter,1); AR_str_all = cell(iter,1);
for k = 1:iter
    BP_str_all{k} = sprintf('[%.3f, %s]', R(k,5), num2str(R(k,6)));
    if isinf(R(k,10)), AR_str_all{k} = 'Inf'; else, AR_str_all{k} = num2str(R(k,10)); end
end
gui_bp_idx = find_bp_idx(BP_OPTS, GT_BP);
gui_match_all = (R(:,3)==GT_POLY & R(:,4)==gui_bp_idx & R(:,7)==GT_SIMULT & ...
                 R(:,8)==GT_MOT24 & R(:,9)==GT_ACOMP & R(:,10)==GT_AROMA & R(:,11)==GT_RMGMR);

T_all = table( (1:iter)', R(:,2), R(:,3), BP_str_all, logical(R(:,7)), logical(R(:,8)), ...
               R(:,9), AR_str_all, logical(R(:,11)), ...
               R(:,12), R(:,13), R(:,14), R(:,15), gui_match_all, S, R(:,16), ...
    'VariableNames', {'Iter','Stage','POLY_ORD','BP_HZ','SIMULT','MOT24', ...
                      'N_ACOMP','N_AROMA','RM_GMR','Validity','Quality', ...
                      'Sensitivity','Mean_QC','GUI_verified','Status','Runtime_s'});
is_new_row = T_all.Stage ~= 0;
T_new = T_all(is_new_row, :);
[~, ord_new] = sort(T_new.Mean_QC, 'descend', 'MissingPlacement','last');
T_new = T_new(ord_new, :);
writetable(T_new, xlsFile, 'Sheet','All Results');

[~, ord_all] = sort(T_all.Mean_QC, 'descend', 'MissingPlacement','last');
T_all_sorted = T_all(ord_all, :);
best_row = T_all_sorted(1, :);
best_row.Origin = string(ternary(best_row.Stage==0,'prior','new'));
writetable(best_row, xlsFile, 'Sheet','Best Pipeline');
writetable(T_all_sorted, xlsFile, 'Sheet','Merged History');
log_msg(logfid, 'Spreadsheet: %s', xlsFile);

rfid = fopen(reportFile, 'w');
fprintf(rfid, '# Denoising Pipeline Optimization Report (with RM_GMR / GSR)\n\n');
fprintf(rfid, '- Generated: %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf(rfid, '- Project:   `%s`\n', projectFile);
fprintf(rfid, '- Subjects:  %d\n', N_SUBJ);
fprintf(rfid, '- Iterations: %d (Stage 1: %d, Stage 2: %d, Prior: %d)\n\n', ...
    iter, sum(R(:,2)==1), sum(R(:,2)==2), n_prior);

fprintf(rfid, '## New Parameter: RM_GMR\n\n');
fprintf(rfid, 'Binary flag. RM_GMR=1 adds a `Grey Matter` entry to `Preproc.confounds` (dim=1, deriv=0, power=1), equivalent to GSR. Motivated by Ciric 2017 / Mascali 2021 benchmarks.\n\n');

fprintf(rfid, '## Validation\n\n');
fprintf(rfid, '| Score | GUI | Re-computed | Diff |\n|---|---|---|---|\n');
fprintf(rfid, '| Validity    | %.4f | %.4f | %+.4f |\n', GT_V,VAL_V,VAL_V-GT_V);
fprintf(rfid, '| Quality     | %.4f | %.4f | %+.4f |\n', GT_Q,VAL_Q,VAL_Q-GT_Q);
fprintf(rfid, '| Sensitivity | %.4f | %.4f | %+.4f |\n\n', GT_S,VAL_S,VAL_S-GT_S);

fprintf(rfid, '## Best Parameter Combination\n\n');
fprintf(rfid, '| Parameter | Value |\n|---|---|\n');
fprintf(rfid, '| POLY_ORD  | %d |\n', best_row.POLY_ORD);
fprintf(rfid, '| BP_HZ     | %s |\n', best_row.BP_HZ{1});
fprintf(rfid, '| SIMULT    | %d |\n', best_row.SIMULT);
fprintf(rfid, '| MOT24     | %d |\n', best_row.MOT24);
fprintf(rfid, '| N_ACOMP   | %d |\n', best_row.N_ACOMP);
fprintf(rfid, '| N_AROMA   | %s |\n', AR_str_all{ord_all(1)});
fprintf(rfid, '| **RM_GMR** | **%d** |\n', best_row.RM_GMR);
fprintf(rfid, '| Validity  | %.4f |\n', best_row.Validity);
fprintf(rfid, '| Quality   | %.4f |\n', best_row.Quality);
fprintf(rfid, '| Sensitivity | %.4f |\n', best_row.Sensitivity);
fprintf(rfid, '| **Mean_QC** | **%.4f** |\n\n', best_row.Mean_QC);

fprintf(rfid, '## Score Distributions (main-effect summary)\n\n');
valid_mask = ~isnan(R(:,15));  Rv = R(valid_mask, :);
param_info = { ...
    'POLY_ORD', 3,  POLY_VALS;
    'BP_HZ',    4,  1:numel(BP_OPTS);
    'SIMULT',   7,  SIMULT_VALS;
    'MOT24',    8,  MOT24_VALS;
    'N_ACOMP',  9,  unique(Rv(:,9))';
    'N_AROMA', 10,  unique(Rv(:,10))';
    'RM_GMR',  11,  RM_GMR_VALS};
for ip = 1:size(param_info,1)
    pname = param_info{ip,1}; pcol = param_info{ip,2}; levels = param_info{ip,3};
    fprintf(rfid, '### %s\n\n', pname);
    fprintf(rfid, '| Level | Mean_QC | Validity | Quality | Sensitivity | n |\n');
    fprintf(rfid, '|---|---|---|---|---|---|\n');
    for lev = levels(:)'
        mask = Rv(:,pcol) == lev;
        if sum(mask)==0, continue; end
        if strcmp(pname,'BP_HZ'), lev_str = mat2str(BP_OPTS{lev});
        elseif isinf(lev),        lev_str = 'Inf';
        else,                     lev_str = num2str(lev);
        end
        fprintf(rfid, '| %s | %.4f | %.4f | %.4f | %.4f | %d |\n', lev_str, ...
            mean(Rv(mask,15),'omitnan'), mean(Rv(mask,12),'omitnan'), ...
            mean(Rv(mask,13),'omitnan'), mean(Rv(mask,14),'omitnan'), sum(mask));
    end
    fprintf(rfid, '\n');
end

fprintf(rfid, '## RM_GMR paired comparison\n\n');
tab_w  = T_all(T_all.RM_GMR==1 & ~isnan(T_all.Mean_QC), :);
tab_wo = T_all(T_all.RM_GMR==0 & ~isnan(T_all.Mean_QC), :);
lifts = [];
for i = 1:height(tab_w)
    r = tab_w(i,:);
    mask = (tab_wo.POLY_ORD==r.POLY_ORD) & strcmp(tab_wo.BP_HZ, r.BP_HZ{1}) & ...
           (tab_wo.SIMULT==r.SIMULT) & (tab_wo.MOT24==r.MOT24) & ...
           (tab_wo.N_ACOMP==r.N_ACOMP) & strcmp(tab_wo.N_AROMA, r.N_AROMA{1});
    if any(mask), lifts(end+1,1) = r.Mean_QC - mean(tab_wo.Mean_QC(mask)); end %#ok<AGROW>
end
if ~isempty(lifts)
    fprintf(rfid, '- Paired cells: %d\n', numel(lifts));
    fprintf(rfid, '- Mean Mean_QC lift (GMR=1 minus GMR=0): %+.4f (std %.4f)\n', mean(lifts), std(lifts));
    fprintf(rfid, '- %% of cells where GMR=1 is better: %.1f%%\n\n', 100*mean(lifts>0));
else
    fprintf(rfid, 'No paired cells available.\n\n');
end

fprintf(rfid, '## Top 10 (RM_GMR=1 arm)\n\n');
T_gmr = T_new(T_new.RM_GMR==1, :);
T_gmr10 = T_gmr(1:min(10,height(T_gmr)), :);
fprintf(rfid, '| POLY | BP | SIMULT | MOT24 | N_ACOMP | N_AROMA | Validity | Quality | Sens | Mean_QC |\n');
fprintf(rfid, '|---|---|---|---|---|---|---|---|---|---|\n');
for r = 1:height(T_gmr10)
    fprintf(rfid, '| %d | %s | %d | %d | %d | %s | %.4f | %.4f | %.4f | %.4f |\n', ...
        T_gmr10.POLY_ORD(r), T_gmr10.BP_HZ{r}, T_gmr10.SIMULT(r), T_gmr10.MOT24(r), ...
        T_gmr10.N_ACOMP(r), T_gmr10.N_AROMA{r}, T_gmr10.Validity(r), T_gmr10.Quality(r), ...
        T_gmr10.Sensitivity(r), T_gmr10.Mean_QC(r));
end

fprintf(rfid, '\n## Top 10 (RM_GMR=0 arm, for comparison)\n\n');
T_nog = T_new(T_new.RM_GMR==0, :);
T_nog10 = T_nog(1:min(10,height(T_nog)), :);
fprintf(rfid, '| POLY | BP | SIMULT | MOT24 | N_ACOMP | N_AROMA | Validity | Quality | Sens | Mean_QC |\n');
fprintf(rfid, '|---|---|---|---|---|---|---|---|---|---|\n');
for r = 1:height(T_nog10)
    fprintf(rfid, '| %d | %s | %d | %d | %d | %s | %.4f | %.4f | %.4f | %.4f |\n', ...
        T_nog10.POLY_ORD(r), T_nog10.BP_HZ{r}, T_nog10.SIMULT(r), T_nog10.MOT24(r), ...
        T_nog10.N_ACOMP(r), T_nog10.N_AROMA{r}, T_nog10.Validity(r), T_nog10.Quality(r), ...
        T_nog10.Sensitivity(r), T_nog10.Mean_QC(r));
end

fprintf(rfid, '\n## Caveats\n\n');
n_fail = sum(cellfun(@(s)startsWith(s,'FAIL'), S));
n_attempted = max(sum(R(:,2) ~= 0), 1);
fprintf(rfid, '- %d of %d new combinations failed (%.1f%%).\n', n_fail, n_attempted, 100*n_fail/n_attempted);
fprintf(rfid, '- Mean_QC = simple mean of Validity, Quality, Sensitivity.\n');
fprintf(rfid, '- RM_GMR implementation: adds Grey Matter confound (dim=1, deriv=0, power=1).\n');
fclose(rfid);
log_msg(logfid, 'Report: %s', reportFile);

fprintf('\n==============================\n');
fprintf(' OPTIMIZATION COMPLETE (RM_GMR)\n');
fprintf('==============================\n');
fprintf(' Iterations : %d  (failed: %d)\n', iter, n_fail);
fprintf(' Best Mean_QC = %.4f\n', best_row.Mean_QC);
fprintf('   POLY=%d  BP=%s  SIM=%d  M24=%d  AC=%d  AR=%s  GMR=%d\n', ...
    best_row.POLY_ORD, best_row.BP_HZ{1}, best_row.SIMULT, best_row.MOT24, ...
    best_row.N_ACOMP, AR_str_all{ord_all(1)}, best_row.RM_GMR);
fprintf(' XLSX   : %s\n', xlsFile);
fprintf(' Report : %s\n', reportFile);
fprintf(' Log    : %s\n', logFile);
fprintf('==============================\n');
log_msg(logfid, '====== OPTIMIZATION FINISHED ======');
fclose(logfid);


%% ====================== LOCAL FUNCTIONS =================================

function [v,q,s,mean_qc,elapsed,status] = run_one_combo( ...
        orig_preproc, poly, bp, simult, mot24, acomp, aroma, rm_gmr, l2_covars, sens_vars)
    global CONN_x;
    v = NaN; q = NaN; s = NaN; mean_qc = NaN; elapsed = 0;
    try
        t0 = tic;
        CONN_x.Preproc = orig_preproc;
        CONN_x.Preproc.detrending = poly;
        CONN_x.Preproc.filter     = bp;
        CONN_x.Preproc.regbp      = 1 + simult;
        CONN_x.Preproc.confounds  = update_confounds(orig_preproc.confounds, mot24, acomp, aroma, rm_gmr);
        v = safe_score(@() conn_qascores('DataValidity', [], []));
        q = safe_score(@() conn_qascores('DataQuality',  [], [], l2_covars, {}));
        s = compute_sensitivity_manual(sens_vars);
        mean_qc = mean([v, q, s], 'omitnan');
        elapsed = toc(t0); status = 'OK';
    catch err
        elapsed = toc(t0);
        status  = ['FAIL: ', strrep(err.message, newline, ' ')];
    end
end

function out = safe_score(fn)
    out = fn(); if isempty(out), out = NaN; end
end

function s = compute_sensitivity_manual(sens_vars)
    global CONN_x;
    folderout = fullfile(CONN_x.folders.qa, ...
        ['QA_GUIrequest_DataSensitivity_', datestr(now,'yyyy_mm_dd')]);
    if ~exist(folderout, 'dir'), mkdir(folderout); end
    conn_process('qaplots', folderout, {'QA_COV'}, 1:CONN_x.Setup.nsubjects, [], [], sens_vars);
    valid = conn_module('get', 'l2covariates', 'QC_ValidSubjects');
    dof   = conn_module('get', 'l2covariates', 'QC_DOF_WelchSatterthwaite');
    if isempty(dof)
        dof = conn_module('get', 'l2covariates', 'QC_DOF');
    end
    m = valid > 0 & ~isnan(valid);
    sum_eff = sum(dof(m) - 3);
    if ~isfinite(sum_eff) || sum_eff < 0
        s = NaN;
    else
        s = 1 - spm_Ncdf(1.645 - 0.1003 * sqrt(sum_eff));
    end
end

function conf = update_confounds(conf_orig, mot24, n_acomp, n_aroma, rm_gmr)
    conf = conf_orig;
    % realignment (24P / 12P)
    idx = find(strcmp(conf.names, 'realignment'), 1);
    if ~isempty(idx)
        if iscell(conf.deriv), conf.deriv{idx} = 1; end
        if iscell(conf.power), conf.power{idx} = 1 + double(logical(mot24)); end
    end
    % aCompCor
    idx = find(strcmp(conf.names, 'aCompCor'), 1);
    if ~isempty(idx) && iscell(conf.dimensions)
        d = conf.dimensions{idx};
        if isempty(d), d = [n_acomp, n_acomp]; else, d(1) = n_acomp; end
        conf.dimensions{idx} = d;
    end
    % aroma
    idx = find(strcmp(conf.names, 'aroma'), 1);
    if n_aroma == 0
        if ~isempty(idx), conf = remove_confound(conf, idx); end
    elseif ~isempty(idx) && iscell(conf.dimensions)
        d = conf.dimensions{idx};
        if isempty(d), d = [n_aroma, n_aroma]; else, d(1) = n_aroma; end
        conf.dimensions{idx} = d;
    elseif isempty(idx)
        conf = append_confound(conf, 'aroma', 'cov', n_aroma);
    end
    % Grey Matter (GSR)  <-- NEW
    idx = find(strcmp(conf.names, 'Grey Matter'), 1);
    if rm_gmr
        if isempty(idx)
            conf = append_confound(conf, 'Grey Matter', 'roi', 1);
        else
            conf.dimensions{idx} = [1, 1];
            if iscell(conf.deriv),  conf.deriv{idx}  = 0; end
            if iscell(conf.power),  conf.power{idx}  = 1; end
            if iscell(conf.filter), conf.filter{idx} = 0; end
        end
    else
        if ~isempty(idx), conf = remove_confound(conf, idx); end
    end
end

function conf = remove_confound(conf, idx)
    flds = fieldnames(conf);
    for f = 1:numel(flds)
        val = conf.(flds{f});
        if iscell(val) && numel(val) >= idx
            val(idx) = [];
            conf.(flds{f}) = val;
        end
    end
end

function conf = append_confound(conf, name, type, dim)
    conf.names{end+1}      = name;
    conf.types{end+1}      = type;
    conf.power{end+1}      = 1;
    conf.deriv{end+1}      = 0;
    conf.dimensions{end+1} = [dim, dim];
    if isfield(conf,'filter'), conf.filter{end+1} = 0; end
    if isfield(conf,'fixed'),  conf.fixed{end+1}  = 0; end
end

function idx = find_bp_idx(bp_opts, target)
    idx = 0;
    for k = 1:numel(bp_opts)
        if isequal(bp_opts{k}, target), idx = k; return; end
    end
end

function idx = find_bp_idx_from_str(bp_opts, s)
    idx = 0;
    if iscell(s), s = s{1}; end
    s = strtrim(char(s));
    for k = 1:numel(bp_opts)
        ref = sprintf('[%.3f, %s]', bp_opts{k}(1), num2str(bp_opts{k}(2)));
        if strcmp(s, ref), idx = k; return; end
    end
end

function v = norm_aroma_val(x)
    if iscell(x), x = x{1}; end
    if isnumeric(x), v = double(x);
    elseif ischar(x) || isstring(x)
        s = strtrim(char(x));
        if strcmpi(s,'inf'), v = Inf; else, v = str2double(s); end
    else, v = NaN; end
end

function k = make_param_key(poly, bp_idx, simult, mot24, acomp, aroma, rm_gmr)
    if isinf(aroma), aroma_str = 'Inf'; else, aroma_str = sprintf('%d', round(aroma)); end
    k = sprintf('%d|%d|%d|%d|%d|%s|%d', round(poly), round(bp_idx), ...
                double(logical(simult)), double(logical(mot24)), ...
                round(acomp), aroma_str, double(logical(rm_gmr)));
end

function [tested_set, prev_T] = load_prior_results(file_list, BP_OPTS, logfid)
    tested_set = containers.Map('KeyType','char','ValueType','logical');
    prev_T = table();
    if isempty(file_list)
        log_msg(logfid, 'EXISTING_RESULT_FILES empty.');
        return;
    end
    for iFile = 1:numel(file_list)
        fpath = file_list{iFile};
        if ~exist(fpath, 'file'), log_msg(logfid, '  WARN not found: %s', fpath); continue; end
        try,    Tprev = readtable(fpath, 'Sheet','All Results');
        catch,  try, Tprev = readtable(fpath);
                catch err, log_msg(logfid,'  WARN %s: %s', fpath, err.message); continue; end
        end
        log_msg(logfid, '  loaded %d rows from %s', height(Tprev), fpath);
        for r = 1:height(Tprev)
            bp_idx = find_bp_idx_from_str(BP_OPTS, Tprev.BP_HZ{r});
            aroma_v = norm_aroma_val(Tprev.N_AROMA(r));
            key = make_param_key(Tprev.POLY_ORD(r), bp_idx, ...
                                 Tprev.SIMULT(r), Tprev.MOT24(r), ...
                                 Tprev.N_ACOMP(r), aroma_v, 0);
            tested_set(key) = true;
        end
        prev_T = vertcat_safe(prev_T, Tprev);
    end
end

function T = vertcat_safe(A, B)
    if isempty(A), T = B; return; end
    if isempty(B), T = A; return; end
    common = intersect(A.Properties.VariableNames, B.Properties.VariableNames, 'stable');
    T = [A(:, common); B(:, common)];
end

function out = ternary(cond, a, b), if cond, out = a; else, out = b; end, end

function log_msg(fid, fmt, varargin)
    msg  = sprintf(fmt, varargin{:});
    ts   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    line = sprintf('[%s] %s\n', ts, msg);
    fprintf('%s', line);
    if fid > 0, fprintf(fid, '%s', line); end
end
