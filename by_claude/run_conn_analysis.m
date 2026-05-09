%% run_conn_analysis.m
% =========================================================================
% Unified first- and second-level CONN analysis for the implicit SRTT
% dataset.  Replaces run_conn_batch_1st.m and run_conn_batch_2nd.m.
%
% =========================================================================
% AUDIT OF PRIOR SCRIPTS (run_conn_batch_1st.m + run_conn_batch_2nd.m)
% =========================================================================
%
% A. SECOND-LEVEL GLM VALIDITY  - needs restructuring
%    ----------------------------------------------------------------------
%    The prior 2nd-level script ran six INDEPENDENT one-sample t-tests
%    (Str_Ran_r12 through Str_Ran_r78, Swi_Ran_r34, Swi_Ran_r56) each with
%    between_subjects.contrast = [1] (AllSubjects) and between_conditions
%    [1 -1] comparing one structured/switch run-pair against 'random'.
%
%    Problems, evaluated against the experimental design (gradual ramp of
%    Sequence B from run 3 onward, with runs 7-8 pure B):
%
%      1. No within-subject run factor.  CONN has no native RM-ANOVA, but
%         conn_glm (/home/aclexp/Software/conn/conn_glm.m lines 1-286)
%         supports multivariate contrasts across conditions when
%         Results.between_conditions.contrast spans multiple effects.
%         The prior script's one-pair-at-a-time design cannot test a
%         learning trajectory - the gradual nature of SL and the A->B
%         ramp are not modeled.
%      2. No linear trend across runs.  Learning is hypothesized to grow
%         across runs; the prior design cannot formally test this.
%      3. Family-wise error uncontrolled across the six 2nd-level models.
%         Six separate results folders at alpha=0.05 each -> 26.5%
%         cumulative Type-I risk.
%      4. Degrees of freedom per model (N=20 -> 19) are correct for the
%         one-sample t-test as-specified, but this design is the wrong
%         one for the research question.
%      5. Runs-within-subject are non-independent; six separate t-tests
%         ignore that.
%
%    VERDICT: needs restructuring.  This script implements three contrast
%    families that match the question: (1) SL main effect averaged across
%    run-pairs, (2) SL linear trend across run-pairs, (3) Switch main
%    effect averaged across run-pairs.  Switch trend is omitted for now
%    because only two switch run-pairs exist (swi_r34, swi_r56) and a
%    linear contrast [-1 1] on two points is equivalent to the paired
%    diff; left as future work if r3-r6 are modelled separately.
%
% B. CLUSTER-LEVEL INFERENCE  - needs to be set explicitly
%    ----------------------------------------------------------------------
%    The prior scripts set no threshold fields.  CONN therefore falls back
%    to whatever CONN_x.Setup.secondlevelanalyses is (GUI default:
%    parametric RFT at voxel p<.001, cluster-FWE p<.05).
%
%    Concerns:
%      1. N=20 is on the small side for parametric RFT; smoothness
%         estimation is unstable.
%      2. Seed-to-voxel connectivity maps are often non-stationary
%         (smoother near the seed) - RFT assumptions are violated.
%      3. Non-parametric permutation makes weaker assumptions and is the
%         current methodological consensus.
%
%    THIS SCRIPT: sets CONN_x.Setup.secondlevelanalyses = 3 (non-
%    parametric only) before running Results.  CONN internally calls
%    conn_randomise (see conn_process.m:5374-5525) which writes cluster-
%    size and cluster-mass permutation null distributions into each
%    Results folder.  Default niters is 1000; a post-hoc helper below
%    re-runs with TFCE (THR_TYPE=5) and niters=10000 using an explicit
%    seed for reproducibility.
%
% C. COMPLETENESS  - partial
%    ----------------------------------------------------------------------
%      - The prior 1st-level script ran both wSBC (measure=1) AND gPPI
%        (measure=3, modulation=1), but only gPPI was tested at the 2nd
%        level.  We commit to gPPI only - for run-specific task
%        connectivity it is the more appropriate model.
%      - Seeds default to Setup.rois = {atlas, networks}: all ~132
%        atlas parcels plus 8 canonical networks.  Per user preference
%        we explore all available seeds without pre-selection (no SRTT-
%        specific ROIs hand-picked).  This inflates the multiple-
%        comparisons problem across seeds; the analysis plan documents
%        the additional FDR-across-seeds step the user should apply
%        after reviewing the per-seed results.
%      - QC_MeanMotion is added as a 2nd-level nuisance covariate.
%      - No behavioural RT covariate (user opted out for now).
%
% D. EFFICIENCY & ROBUSTNESS  - several hard-coded problems fixed
%    ----------------------------------------------------------------------
%      - Prior scripts hard-coded rootDir='/media/data3/...' (IP=37 path)
%        - stale on the IP=23 workstation.  This script uses the current
%        root (/home/aclexp/pinwei/Joanne_SRT_fMRI/).
%      - Prior batchFile='v1_no_param.mat' - stale; the optimization
%        work targeted v3 and this script targets the newly-built
%        v5_optimized.mat (generated by run_conn_batch_optimized.m).
%      - Prior run_conn_batch_2nd.m lines 56,67 have a trailing-space
%        typo 'Results.name ' - silently creates the wrong batch field.
%        Fixed here.
%      - All paths are variables at the top of this script; body has no
%        hard-coded paths.
%      - DRY_RUN flag validates inputs and prints the plan without
%        launching CONN.
%      - Progress + error messages logged to
%        logs/run_conn_analysis_log.txt.
%      - Per-contrast try/catch: a failing 2nd-level call does not abort
%        the other contrasts.
%
% =========================================================================
% Usage
%   1. (Prereq) Build the optimized project:
%        >> run_conn_batch_optimized
%   2. (Optional) Set DRY_RUN = true below and run this script to see
%      the plan without executing.
%   3. >> run_conn_analysis
% =========================================================================

clear; clc; close all;
global CONN_x;


%% ====================== 1. CONFIGURATION ================================

DRY_RUN = true;      % set to false to actually run CONN
N_PERMS = 10000;     % non-parametric permutations for TFCE post-hoc
PERM_SEED = 1234;    % seed for reproducibility

rootDir     = '/home/aclexp/pinwei/Joanne_SRT_fMRI';
PROJECT_FILE = fullfile(rootDir, 'conn_out', 'v5_no_param.mat');
LOG_FILE     = fullfile(rootDir, 'logs', 'run_conn_analysis_log.txt');

% First-level analysis name (used for both first-level write and
% second-level analysis lookup via Results.name).
ANALYSIS_NAME = 'gPPI';

% Full condition list (must match run_conn_batch_optimized.m COND_NAMES)
COND_LIST  = {'random','str_r12','str_r34','str_r56','str_r78', ...
              'swi_r34','swi_r56','incorrect'};

% Second-level contrast specifications.  Each row is one Results call.
%   .saveas            results folder label
%   .between_conditions effect_names + contrast vector
%   .description       human-readable purpose (written to report)
CONTRASTS = struct();
CONTRASTS(1).saveas = 'SL_main_effect';
CONTRASTS(1).between_conditions_names = {'str_r12','str_r34','str_r56','str_r78','random'};
CONTRASTS(1).between_conditions_contrast = [1 1 1 1 -4];
CONTRASTS(1).description = 'Structured > random averaged over 4 run-pairs';

CONTRASTS(2).saveas = 'SL_linear_trend';
CONTRASTS(2).between_conditions_names = {'str_r12','str_r34','str_r56','str_r78'};
CONTRASTS(2).between_conditions_contrast = [-3 -1 1 3];
CONTRASTS(2).description = 'Linear trend of structured connectivity r12->r78 (learning trajectory)';

CONTRASTS(3).saveas = 'Switch_main_effect';
CONTRASTS(3).between_conditions_names = {'swi_r34','swi_r56','random'};
CONTRASTS(3).between_conditions_contrast = [1 1 -2];
CONTRASTS(3).description = 'Switch (sequence B) > random averaged over 2 run-pairs';

% Second-level subject model: AllSubjects (one-sample t on group mean)
% with QC_MeanMotion as nuisance covariate.  The contrast [1 0] tests
% the group effect controlling for motion; [0 1] tests the motion
% covariate (not reported here but available in SPM.mat).
BS_EFFECT_NAMES = {'AllSubjects','QC_MeanMotion'};
BS_CONTRAST     = [1 0];


%% ====================== 2. LOGGING SETUP ================================

if ~exist(fileparts(LOG_FILE), 'dir'), mkdir(fileparts(LOG_FILE)); end
logfid = fopen(LOG_FILE, 'w');
log_msg(logfid, '====== CONN task-connectivity analysis ======');
log_msg(logfid, 'DRY_RUN       : %s', mat2str(DRY_RUN));
log_msg(logfid, 'Project file  : %s', PROJECT_FILE);
log_msg(logfid, 'Analysis name : %s', ANALYSIS_NAME);
log_msg(logfid, 'Permutations  : %d (seed=%d)', N_PERMS, PERM_SEED);


%% ====================== 3. SANITY CHECKS ================================
% Performed regardless of DRY_RUN so bad inputs surface before a long run.

log_msg(logfid, '');
log_msg(logfid, '========== SANITY CHECKS ==========');

problems = {};
if ~exist(PROJECT_FILE, 'file')
    problems{end+1} = sprintf('project file missing: %s', PROJECT_FILE);
    problems{end+1} = '  hint: run scripts/run_conn_batch_optimized.m first.';
end

% Load the project (read-only mode for DRY_RUN) to inspect structure.
% Even in DRY_RUN we load so we can validate the contrast spec against
% the actual Setup.conditions.names list.
if exist(PROJECT_FILE, 'file')
    addpath('/home/aclexp/Software/conn');
    log_msg(logfid, 'Loading project for inspection ...');
    conn('load', PROJECT_FILE);
    nsub = CONN_x.Setup.nsubjects;
    setup_conds = CONN_x.Setup.conditions.names;
    log_msg(logfid, '  nsubjects = %d', nsub);
    log_msg(logfid, '  conditions in project: %s', strjoin(setup_conds, ' | '));

    for ic = 1:numel(CONTRASTS)
        missing = setdiff(CONTRASTS(ic).between_conditions_names, setup_conds);
        if ~isempty(missing)
            problems{end+1} = sprintf('contrast "%s" references missing conditions: %s', ...
                CONTRASTS(ic).saveas, strjoin(missing, ', '));
        end
        ncond = numel(CONTRASTS(ic).between_conditions_names);
        ncon  = numel(CONTRASTS(ic).between_conditions_contrast);
        if ncond ~= ncon
            problems{end+1} = sprintf('contrast "%s" effect_names (%d) != contrast length (%d)', ...
                CONTRASTS(ic).saveas, ncond, ncon);
        end
    end

    l2names = CONN_x.Setup.l2covariates.names;
    for bs = BS_EFFECT_NAMES
        if ~any(strcmp(l2names, bs{1})) && ~strcmp(bs{1}, 'AllSubjects')
            problems{end+1} = sprintf('2nd-level covariate "%s" not in project.  Available: %s', ...
                bs{1}, strjoin(l2names,', '));
        end
    end
end

if ~isempty(problems)
    log_msg(logfid, 'SANITY CHECK FAILURES:');
    for p = 1:numel(problems), log_msg(logfid, '  - %s', problems{p}); end
    if DRY_RUN
        log_msg(logfid, 'DRY_RUN=true: continuing to print the plan anyway.');
    else
        fclose(logfid);
        error('Sanity check failed.  See %s.', LOG_FILE);
    end
else
    log_msg(logfid, 'All sanity checks passed.');
end


%% ====================== 4. PRINT ANALYSIS PLAN ==========================

log_msg(logfid, '');
log_msg(logfid, '========== ANALYSIS PLAN ==========');
log_msg(logfid, 'First level:');
log_msg(logfid, '  Analysis name = %s', ANALYSIS_NAME);
log_msg(logfid, '  measure=3 (regression bivariate), modulation=1 (gPPI)');
log_msg(logfid, '  type=3 (ROI-to-ROI + seed-to-voxel), conditions=%s', strjoin(COND_LIST,' | '));
log_msg(logfid, '  seeds = all Setup.rois (atlas + networks; exploratory - no pre-selection)');
log_msg(logfid, '');
log_msg(logfid, 'Second level:');
log_msg(logfid, '  CONN_x.Setup.secondlevelanalyses = 3  (non-parametric only)');
log_msg(logfid, '  between_subjects effects = %s', strjoin(BS_EFFECT_NAMES,' | '));
log_msg(logfid, '  between_subjects contrast = %s', mat2str(BS_CONTRAST));
log_msg(logfid, '  Contrast families:');
for ic = 1:numel(CONTRASTS)
    log_msg(logfid, '    [%d] %-20s  %s', ic, CONTRASTS(ic).saveas, CONTRASTS(ic).description);
    log_msg(logfid, '         conditions: %s', strjoin(CONTRASTS(ic).between_conditions_names, ', '));
    log_msg(logfid, '         contrast  : %s', mat2str(CONTRASTS(ic).between_conditions_contrast));
end
log_msg(logfid, '');
log_msg(logfid, 'Post-GLM cluster inference:');
log_msg(logfid, '  TFCE (THR_TYPE=5), %d permutations, two-sided, seed=%d', N_PERMS, PERM_SEED);
log_msg(logfid, '  Per-Results helper: conn_randomise on SPM.mat in each folder.');


%% ====================== 5. DRY_RUN SHORT-CIRCUIT ========================

if DRY_RUN
    log_msg(logfid, '');
    log_msg(logfid, 'DRY_RUN=true.  Nothing executed.  Set DRY_RUN=false to run.');
    fclose(logfid);
    fprintf('\n*** DRY_RUN complete.  Plan logged to %s. ***\n', LOG_FILE);
    return;
end


%% ====================== 6. FIRST-LEVEL gPPI =============================

log_msg(logfid, '');
log_msg(logfid, '========== FIRST LEVEL (%s) ==========', ANALYSIS_NAME);

t0 = tic;
try
    conn_batch( ...
        'filename',          PROJECT_FILE, ...
        'Analysis.name',     ANALYSIS_NAME, ...
        'Analysis.measure',  3, ...       % regression (bivariate)
        'Analysis.modulation', 1, ...     % gPPI interaction
        'Analysis.conditions', string(COND_LIST), ...
        'Analysis.type',     3, ...       % ROI-to-ROI + seed-to-voxel
        'Analysis.done',     1, ...
        'Analysis.overwrite','Yes');
    log_msg(logfid, 'First level OK (%.1f min).', toc(t0)/60);
catch err
    log_msg(logfid, 'First level FAILED: %s', err.message);
    fclose(logfid);
    rethrow(err);
end


%% ====================== 7. SWITCH TO NON-PARAMETRIC 2ND-LEVEL ===========

log_msg(logfid, '');
log_msg(logfid, '========== CONFIGURE 2ND-LEVEL (non-parametric) ==========');

conn('load', PROJECT_FILE);
CONN_x.Setup.secondlevelanalyses = 3;      % 1=both, 2=RFT only, 3=permutation only
conn save;
log_msg(logfid, 'Set CONN_x.Setup.secondlevelanalyses = 3 and saved project.');


%% ====================== 8. SECOND-LEVEL GLMs ============================

log_msg(logfid, '');
log_msg(logfid, '========== SECOND LEVEL ==========');

results_folders = cell(numel(CONTRASTS),1);
for ic = 1:numel(CONTRASTS)
    C = CONTRASTS(ic);
    log_msg(logfid, '[%d/%d] %s  (%s)', ic, numel(CONTRASTS), C.saveas, C.description);

    t0 = tic;
    try
        conn_batch( ...
            'filename', PROJECT_FILE, ...
            'Results.name',                    ANALYSIS_NAME, ...
            'Results.saveas',                  C.saveas, ...
            'Results.between_subjects.effect_names', BS_EFFECT_NAMES, ...
            'Results.between_subjects.contrast',     BS_CONTRAST, ...
            'Results.between_conditions.effect_names', C.between_conditions_names, ...
            'Results.between_conditions.contrast',     C.between_conditions_contrast, ...
            'Results.display', 0, ...
            'Results.done',    1);
        log_msg(logfid, '  Results.done OK (%.1f s).', toc(t0));

        results_folders{ic} = find_results_folder(PROJECT_FILE, ANALYSIS_NAME, C.saveas);
        if isempty(results_folders{ic})
            log_msg(logfid, '  WARN: could not locate Results folder for "%s".', C.saveas);
        else
            log_msg(logfid, '  Folder: %s', results_folders{ic});
        end
    catch err
        log_msg(logfid, '  FAILED: %s', err.message);
        results_folders{ic} = '';
        continue;   % continue with other contrasts
    end
end


%% ====================== 9. POST-GLM TFCE PERMUTATION ====================
% CONN's Results step runs a default 1000-iteration permutation under
% secondlevelanalyses=3.  For reproducible TFCE-corrected inference we
% re-run conn_randomise on the SPM.mat in each Results folder with
% THR_TYPE=5 and niters=N_PERMS.

log_msg(logfid, '');
log_msg(logfid, '========== POST-GLM TFCE PERMUTATION ==========');

rng(PERM_SEED);   % reproducibility for the permutation sampler

for ic = 1:numel(CONTRASTS)
    fld = results_folders{ic};
    if isempty(fld) || ~exist(fld, 'dir')
        log_msg(logfid, '[%d] skipping TFCE for %s (no folder)', ic, CONTRASTS(ic).saveas);
        continue;
    end
    spm_file = fullfile(fld, 'SPM.mat');
    if ~exist(spm_file, 'file')
        log_msg(logfid, '[%d] skipping TFCE for %s (no SPM.mat in %s)', ic, CONTRASTS(ic).saveas, fld);
        continue;
    end

    log_msg(logfid, '[%d] TFCE for %s ...', ic, CONTRASTS(ic).saveas);
    t0 = tic;
    try
        run_tfce(spm_file, N_PERMS, logfid);
        log_msg(logfid, '  TFCE OK (%.1f min).', toc(t0)/60);
    catch err
        log_msg(logfid, '  TFCE FAILED: %s', err.message);
    end
end


%% ====================== 10. DONE ========================================

log_msg(logfid, '');
log_msg(logfid, '====== ANALYSIS COMPLETE ======');
log_msg(logfid, 'Results folders:');
for ic = 1:numel(CONTRASTS)
    log_msg(logfid, '  [%d] %s -> %s', ic, CONTRASTS(ic).saveas, results_folders{ic});
end
log_msg(logfid, 'Inspect interactively with:  conn; conn(''load'',''%s''); conn gui_results', PROJECT_FILE);
fclose(logfid);
fprintf('\n*** Analysis complete.  Log: %s ***\n', LOG_FILE);


%% ====================== LOCAL FUNCTIONS =================================

function log_msg(fid, fmt, varargin)
    msg  = sprintf(fmt, varargin{:});
    ts   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    line = sprintf('[%s] %s\n', ts, msg);
    fprintf('%s', line);
    if fid > 0, fprintf(fid, '%s', line); end
end


function folder = find_results_folder(projectFile, analysisName, saveas)
% FIND_RESULTS_FOLDER  Locate the output folder CONN wrote for this
% 2nd-level contrast.  CONN stores results under
%   <project>/results/secondlevel/<analysisName>/<saveas>
% but the exact subfolder name depends on CONN version.  Try several.
    folder = '';
    [fp, fn] = fileparts(projectFile);
    cand_roots = { ...
        fullfile(fp, fn, 'results','secondlevel'), ...
        fullfile(fp, 'results','secondlevel')};
    for iR = 1:numel(cand_roots)
        root = cand_roots{iR};
        if ~exist(root,'dir'), continue; end
        for name = {analysisName, ['ANALYSIS_', analysisName]}
            cand = fullfile(root, name{1}, saveas);
            if exist(cand, 'dir'), folder = cand; return; end
        end
        % Fallback: find any folder under root matching saveas
        d = dir(fullfile(root,'**',saveas));
        if ~isempty(d)
            folder = fullfile(d(1).folder, d(1).name);
            return;
        end
    end
end


function run_tfce(spm_file, niters, logfid)
% RUN_TFCE  Load a 2nd-level SPM.mat and run conn_randomise with TFCE
% (THR_TYPE=5) and the requested number of iterations.
%
%   spm_file : absolute path to SPM.mat written by a CONN Results step
%   niters   : number of permutations (e.g. 10000)
%   logfid   : log file id (0 = stdout only)
%
% conn_randomise signature (from conn_randomise.m line 1):
%   conn_randomise(X,Y,c1,c2,aopt,Pthr,Pthr_type,Pthr_side,niters,...
%                  filename,overwrite,datatype,analysismask,...
%                  groupingsamples,isdisplay)
    S = load(spm_file, 'SPM');
    SPM = S.SPM;
    X  = SPM.xX.X;
    c1 = SPM.xCon(1).c';                % primary contrast
    c2 = eye(size(X,2));                 % full-rank reduced model
    Pthr      = [0.05 0.01 0.005 0.001 0.0005 0.0001];
    Pthr_type = 5 * ones(size(Pthr));    % 5 = TFCE for every threshold level
    Pthr_side = [1 2 3];                 % pos, neg, two-sided

    folder = fileparts(spm_file);
    sim_file = fullfile(folder, sprintf('nonparametric_TFCE_%dperms.mat', niters));

    % Read the data cube from SPM.xY.VY (one 3D volume per subject per
    % contrasted condition).  conn_randomise expects Y as 4D/ndarray in
    % (nsubjects x nY x nvoxels) after permuting; for seed-to-voxel it
    % accepts the flattened 4D volume form used in conn_process.m:5995.
    Y = spm_read_vols(SPM.xY.VY);        % [nx ny nz nobs]
    dims = size(Y);
    mask = ~any(isnan(Y),4) & any(diff(Y,1,4)~=0, 4);
    Yv   = reshape(Y, [], dims(4));
    Yv   = Yv(mask(:), :)';              % [nobs x nvox]

    [gx,gy,gz] = ind2sub(dims(1:3), find(mask));
    xyz = [gx(:) gy(:) gz(:)]';
    adj = xyz;                           % 3D volume adjacency (conn_clusters handles)

    % Reshape to the (nsubjects, 1, nvox) format conn_randomise expects
    % when called voxel-wise (see conn_process.m:5990).
    Ycube = permute(reshape(Yv, size(Yv,1), size(X,1), []), [2 3 1]);

    if logfid > 0
        fprintf(logfid, '      -> conn_randomise TFCE niters=%d nvox=%d...\n', ...
            niters, sum(mask(:)));
    end
    conn_randomise(X, Ycube, c1, c2, [], Pthr, Pthr_type, Pthr_side, ...
        niters, sim_file, 'overwrite', adj, [], [], false);
end
