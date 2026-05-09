% fix_l2covariates.m
%
% Surgically restores second-level (l2) covariates in the CONN project .mat
% after they were wiped by an accidental `conn_batch` call that invoked the
% Setup step with `Setup.subjects.groups_names` and `Setup.done=1`.
%
% The script does NOT touch Setup, Denoising, ROIs, conditions, or
% first-level analyses. It only rewrites CONN_x.Setup.l2covariates.
%
% Restored entries (appended after any surviving tissue-volume QC):
%   QC_MeanMotion            - recomputed from fmriprep fd_run-*.tsv
%   QC_InvalidScans          - recomputed from fmriprep scrubbing_run-*.tsv
%   AllSubjects              - 1 for every subject
%   ExcludeOutlierSubjects   - 0 for every subject (no outliers flagged)
%   ' '                      - trailing placeholder required by CONN
%
% A timestamped backup copy of the project .mat is written to conn_out/
% before the in-place rewrite.
%
% Usage (from MATLAB, with no CONN project currently loaded):
%   run('/home/aclexp/pinwei/Joanne_SRT_fMRI/scripts/fix_l2covariates.m')

clear; clc;

%% Paths

rootDir   = '/home/aclexp/pinwei/Joanne_SRT_fMRI/';
connDir   = '/home/aclexp/Software/conn';
spmDir    = '/home/aclexp/Software/spm';
addpath(connDir);
addpath(spmDir);

batchFile = fullfile(rootDir, 'conn_out', 'v5_no_param.mat');
bidsDir   = fullfile(rootDir, 'data', 'fmriprep');

assert(exist(batchFile, 'file') == 2, 'Project .mat not found: %s', batchFile);

%% Subject + session configuration (mirrors run_conn_preproc.m)

SID_LIST = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25];
N_SUBJS  = numel(SID_LIST);
N_RUNS   = 8;

%% Load project and snapshot current state

fprintf('[fix] Loading %s\n', batchFile);
S      = load(batchFile);
assert(isfield(S, 'CONN_x'), 'CONN_x missing from .mat');
CONN_x = S.CONN_x;

assert(CONN_x.Setup.nsubjects == N_SUBJS, ...
    'nsubjects mismatch: project has %d, script expects %d', ...
    CONN_x.Setup.nsubjects, N_SUBJS);

existingNames = CONN_x.Setup.l2covariates.names;
fprintf('[fix] Existing l2covariates (%d):\n', numel(existingNames));
for i = 1:numel(existingNames)
    fprintf('    [%d] %s\n', i, existingNames{i});
end

%% Timestamped backup

ts         = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
backupFile = fullfile(rootDir, 'conn_out', ...
    sprintf('v5_no_param_pre_fix_%s.mat', ts));
copyfile(batchFile, backupFile);
fprintf('[fix] Backup written: %s\n', backupFile);

%% Recompute QC_MeanMotion and QC_InvalidScans from fmriprep covariate tsvs

meanMotion   = zeros(N_SUBJS, N_RUNS);
invalidScans = zeros(N_SUBJS, N_RUNS);

for isub = 1:N_SUBJS
    sid  = SID_LIST(isub);
    subj = sprintf('sub-%03d', sid);
    for irun = 1:N_RUNS
        fdFile    = fullfile(bidsDir, subj, 'func', 'covariates', ...
            sprintf('fd_run-%02d.tsv', irun));
        scrubFile = fullfile(bidsDir, subj, 'func', 'covariates', ...
            sprintf('scrubbing_run-%02d.tsv', irun));
        assert(exist(fdFile, 'file') == 2,    'Missing %s', fdFile);
        assert(exist(scrubFile, 'file') == 2, 'Missing %s', scrubFile);
        Tfd    = readtable(fdFile,    'FileType', 'text', 'Delimiter', '\t');
        Tscrub = readtable(scrubFile, 'FileType', 'text', 'Delimiter', '\t');
        meanMotion(isub, irun)   = mean(Tfd{:, 1});
        invalidScans(isub, irun) = sum(Tscrub{:, 1});
    end
end

QC_MeanMotion   = mean(meanMotion, 2);
QC_InvalidScans = sum(invalidScans, 2);

fprintf('[fix] QC_MeanMotion   : min=%.4f max=%.4f mean=%.4f\n', ...
    min(QC_MeanMotion), max(QC_MeanMotion), mean(QC_MeanMotion));
fprintf('[fix] QC_InvalidScans : min=%d max=%d mean=%.1f\n', ...
    min(QC_InvalidScans), max(QC_InvalidScans), mean(QC_InvalidScans));

%% Build new l2covariates list
%
% Keep every surviving non-empty entry (tissue volumes), drop the trailing
% empty ' ' placeholder, then append the four restored entries and the
% placeholder at the end.

baseNames   = existingNames;
baseDescrip = CONN_x.Setup.l2covariates.descrip;
if numel(baseDescrip) < numel(baseNames)
    baseDescrip(end+1:numel(baseNames)) = {''};
elseif numel(baseDescrip) > numel(baseNames)
    baseDescrip = baseDescrip(1:numel(baseNames));
end

isPlaceholder = cellfun(@(x) isempty(strtrim(x)), baseNames);
baseNames   = baseNames(~isPlaceholder);
baseDescrip = baseDescrip(~isPlaceholder);

addNames = {'QC_MeanMotion', 'QC_InvalidScans', 'AllSubjects', 'ExcludeOutlierSubjects'};
addDescrip = repmat({''}, 1, numel(addNames));

newNames   = [baseNames,   addNames,   {' '}];
newDescrip = [baseDescrip, addDescrip, {''}];

nBase = numel(baseNames);
nAll  = numel(newNames);

%% Build per-subject values

newValues = cell(1, N_SUBJS);
for isub = 1:N_SUBJS
    oldVals = {};
    if isub <= numel(CONN_x.Setup.l2covariates.values)
        oldVals = CONN_x.Setup.l2covariates.values{isub};
    end

    vals = cell(1, nAll);
    % Preserve base (tissue-volume) values in their original slots.
    baseIdx = find(~isPlaceholder);
    for k = 1:numel(baseIdx)
        src = baseIdx(k);
        if src <= numel(oldVals), vals{k} = oldVals{src}; else, vals{k} = []; end
    end
    % Append restored entries.
    vals{nBase + 1} = QC_MeanMotion(isub);
    vals{nBase + 2} = QC_InvalidScans(isub);
    vals{nBase + 3} = 1;   % AllSubjects
    vals{nBase + 4} = 0;   % ExcludeOutlierSubjects (no subject flagged)
    vals{nBase + 5} = [];  % trailing placeholder

    newValues{isub} = vals;
end

CONN_x.Setup.l2covariates.names   = newNames;
CONN_x.Setup.l2covariates.descrip = newDescrip;
CONN_x.Setup.l2covariates.values  = newValues;

%% Save back (match original v7 format)

save(batchFile, 'CONN_x', '-v7');
fprintf('[fix] Patched project saved: %s\n', batchFile);

%% Verification summary

fprintf('\n[fix] Resulting l2covariates (%d):\n', numel(newNames));
for i = 1:numel(newNames)
    if i <= nBase
        label = '(preserved)';
    elseif i == nAll
        label = '(placeholder)';
    else
        label = '(restored)';
    end
    vals_i = nan(N_SUBJS, 1);
    for isub = 1:N_SUBJS
        v = newValues{isub}{i};
        if ~isempty(v) && isnumeric(v), vals_i(isub) = v(1); end
    end
    fprintf('    [%2d] %-28s %-13s min=%-10.4g max=%-10.4g mean=%-10.4g\n', ...
        i, newNames{i}, label, min(vals_i), max(vals_i), mean(vals_i, 'omitnan'));
end

fprintf('\n[fix] Done. AllSubjects and ExcludeOutlierSubjects are now present.\n');
