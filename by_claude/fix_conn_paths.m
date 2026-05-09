% fix_conn_paths.m
%
% Rebase absolute paths inside a CONN project .mat that was copied from
% another machine. Walks every string in the file (CONN_x struct + any
% side-car .mat files under the project folder) and applies one or more
% prefix replacements, e.g. /media/data3/Joanne_SRT_pw/ -> /home/aclexp/pinwei/Joanne_SRT_fMRI/.
%
% USAGE
%   1. Confirm RSYNC (or whatever copied the project) has FINISHED before running.
%   2. Keep DRY_RUN = true and run once. Inspect the log — confirm the old roots
%      and planned replacements.
%   3. If everything looks right, set DRY_RUN = false and run again.
%   4. Open the project in CONN to verify files resolve.

%% -------- CONFIG ----------------------------------------------------------
PROJECT_MAT    = '/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/no_PM_260421.mat';
PROJECT_FOLDER = '/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/no_PM_260421';

% ordered replacement list: {OLD_PREFIX, NEW_PREFIX}
% First match wins; be specific — put longer/more-specific prefixes first.
REPLACEMENTS = {
    '/media/data3/Joanne_SRT_pw', '/home/aclexp/pinwei/Joanne_SRT_fMRI';
};

LOG_PATH = '/home/aclexp/pinwei/Joanne_SRT_fMRI/temp/fix_conn_paths.log';
DRY_RUN  = false;  % set false only after reviewing the dry-run log
MAX_PREVIEW_PER_FILE = 8;  % how many (before/after) pairs to print per .mat
%% --------------------------------------------------------------------------

if ~exist(fileparts(LOG_PATH),'dir'), mkdir(fileparts(LOG_PATH)); end
fid = fopen(LOG_PATH,'w');
logf = @(varargin) multi_fprintf(fid, varargin{:});

logf('fix_conn_paths - %s\n', datestr(now));
logf('DRY_RUN = %d\n', DRY_RUN);
logf('Replacements:\n');
for r = 1:size(REPLACEMENTS,1)
    logf('  %-55s  ->  %s\n', REPLACEMENTS{r,1}, REPLACEMENTS{r,2});
end
logf('\n');

% Build list of .mat files to process
files = {PROJECT_MAT};
if isfolder(PROJECT_FOLDER)
    d = dir(fullfile(PROJECT_FOLDER,'**','*.mat'));
    for i = 1:numel(d)
        p = fullfile(d(i).folder, d(i).name);
        if contains(p, '.backup'), continue; end
        files{end+1,1} = p; %#ok<SAGROW>
    end
end
logf('Files to process: %d\n', numel(files));

totals = struct('rewrote',0,'unmatched_abs',0);
all_unmatched = {};
prefix_hist = containers.Map('KeyType','char','ValueType','int64');

for i = 1:numel(files)
    [r, u, h] = process_one(files{i}, REPLACEMENTS, DRY_RUN, MAX_PREVIEW_PER_FILE, logf);
    totals.rewrote       = totals.rewrote       + r;
    totals.unmatched_abs = totals.unmatched_abs + numel(u);
    all_unmatched        = [all_unmatched; u(:)];
    ks = keys(h);
    for k = 1:numel(ks)
        key = ks{k};
        if isKey(prefix_hist, key)
            prefix_hist(key) = prefix_hist(key) + h(key);
        else
            prefix_hist(key) = h(key);
        end
    end
end

logf('\n===== Summary =====\n');
logf('Files processed       : %d\n', numel(files));
logf('Strings rewritten     : %d\n', totals.rewrote);
logf('Unmatched abs paths   : %d unique\n', numel(unique(all_unmatched)));

logf('\n===== Absolute-path prefix histogram (first 4 components) =====\n');
ks = keys(prefix_hist);
counts = zeros(1, numel(ks));
for k = 1:numel(ks), counts(k) = prefix_hist(ks{k}); end
[~, ord] = sort(counts,'descend');
for k = 1:numel(ord)
    key = ks{ord(k)};
    logf('  %8d   %s\n', prefix_hist(key), key);
end

logf('\n===== Unmatched absolute paths (sample, up to 50) =====\n');
uu = unique(all_unmatched);
for i = 1:min(numel(uu), 50), logf('  %s\n', uu{i}); end
if numel(uu) > 50, logf('  ... (%d more)\n', numel(uu)-50); end

fclose(fid);
fprintf('\nDone. DRY_RUN=%d. Log: %s\n', DRY_RUN, LOG_PATH);


%% ========================================================================
function [n_rewrote, unmatched_list, prefix_hist] = process_one(mat_path, reps, dry_run, preview_cap, logf)
    logf('\n--- %s\n', mat_path);
    n_rewrote = 0; unmatched_list = {};
    prefix_hist = containers.Map('KeyType','char','ValueType','int64');

    if ~isfile(mat_path)
        logf('  [MISSING]\n'); return;
    end

    backup = [mat_path '.backup'];
    if ~dry_run && ~isfile(backup)
        copyfile(mat_path, backup);
        logf('  backup -> %s\n', backup);
    end

    try
        S = load(mat_path);
    catch ME
        logf('  [LOAD FAILED] %s\n', ME.message);
        return;
    end

    state = struct('reps',{reps},'n',0,'unmatched',{{}}, ...
                   'hist',prefix_hist,'preview',{{}},'preview_cap',preview_cap);
    fns = fieldnames(S);
    for i = 1:numel(fns)
        [S.(fns{i}), state] = walk(S.(fns{i}), state);
    end

    for k = 1:numel(state.preview)
        logf('    %s\n', state.preview{k});
    end

    n_rewrote = state.n;
    unmatched_list = state.unmatched;
    prefix_hist = state.hist;
    logf('  rewrote %d, unmatched abs %d\n', n_rewrote, numel(unmatched_list));

    if ~dry_run && n_rewrote > 0
        v = matfile_version(mat_path);
        save(mat_path, '-struct', 'S', v);
        logf('  saved (%s)\n', v);
    elseif dry_run
        logf('  [DRY_RUN] not saved\n');
    else
        logf('  no changes; not saved\n');
    end
end

% -------------------------------------------------------------------------
function [x, state] = walk(x, state)
    if iscell(x)
        for i = 1:numel(x)
            [x{i}, state] = walk(x{i}, state);
        end
    elseif isstruct(x)
        if numel(x) > 1
            for i = 1:numel(x)
                xi = x(i);
                [xi, state] = walk(xi, state);
                x(i) = xi;
            end
        elseif numel(x) == 1
            fn = fieldnames(x);
            for k = 1:numel(fn)
                [x.(fn{k}), state] = walk(x.(fn{k}), state);
            end
        end
    elseif ischar(x) && ~isempty(x)
        [x, state] = rebase_char(x, state);
    elseif isstring(x)
        for i = 1:numel(x)
            [s2, state] = rebase_char(char(x(i)), state);
            x(i) = string(s2);
        end
    end
end

% -------------------------------------------------------------------------
function [out, state] = rebase_char(s, state)
    lines = regexp(s, '\r?\n', 'split');
    any_changed = false;
    for i = 1:numel(lines)
        ln = lines{i};
        if ~looks_like_path(ln), continue; end
        if is_abs(ln)
            key = prefix4(ln);
            if isKey(state.hist, key)
                state.hist(key) = state.hist(key) + 1;
            else
                state.hist(key) = 1;
            end
        end
        new_ln = apply_reps(ln, state.reps);
        if ~strcmp(ln, new_ln)
            lines{i} = new_ln;
            any_changed = true;
            state.n = state.n + 1;
            if numel(state.preview) < state.preview_cap
                state.preview{end+1} = sprintf('- %s\n    + %s', ln, new_ln);
            end
        elseif is_abs(ln) && ~starts_with_any(ln, state.reps(:,2))
            state.unmatched{end+1,1} = ln;
        end
    end
    if any_changed
        out = strjoin(lines, char(10));
    else
        out = s;
    end
end

% -------------------------------------------------------------------------
function new_s = apply_reps(s, reps)
    new_s = s;
    for r = 1:size(reps,1)
        old_p = reps{r,1}; new_p = reps{r,2};
        if endsWith(old_p,'/') || endsWith(old_p,'\')
            old_p = old_p(1:end-1);
        end
        if endsWith(new_p,'/') || endsWith(new_p,'\')
            new_p = new_p(1:end-1);
        end
        if startsWith(new_s, [old_p '/']) || startsWith(new_s, [old_p '\']) || strcmp(new_s, old_p)
            new_s = [new_p new_s(numel(old_p)+1:end)];
            new_s = strrep(new_s, '\', '/');
            return;
        end
    end
end

function tf = starts_with_any(s, prefs)
    tf = false;
    for i = 1:numel(prefs)
        p = prefs{i};
        if endsWith(p,'/') || endsWith(p,'\'), p = p(1:end-1); end
        if startsWith(s, [p '/']) || strcmp(s, p)
            tf = true; return;
        end
    end
end

function tf = looks_like_path(s)
    tf = (contains(s,'/') || contains(s,'\')) && ( is_abs(s) || ...
         contains(s,'.nii') || contains(s,'.mat') || contains(s,'.img') || ...
         contains(s,'.txt') || contains(s,'.tsv') || contains(s,'.json') || ...
         contains(s,'.gz')  || contains(s,'.hdr') || contains(s,'.nii.gz') );
end

function tf = is_abs(s)
    tf = false;
    if isempty(s), return; end
    c1 = s(1);
    if c1 == '/' || c1 == '\'
        tf = true; return;
    end
    if numel(s) >= 3 && isletter(s(1)) && s(2) == ':' && (s(3) == '\' || s(3) == '/')
        tf = true;
    end
end

function p = prefix4(s)
    parts = regexp(s, '[/\\]', 'split');
    parts = parts(~cellfun('isempty', parts));
    n = min(4, numel(parts));
    p = ['/' strjoin(parts(1:n), '/')];
end

function v = matfile_version(p)
    fid = fopen(p,'r'); header = fread(fid, 128, '*char')'; fclose(fid);
    if startsWith(header, 'MATLAB 7.3')
        v = '-v7.3';
    elseif startsWith(header, 'MATLAB 5.0 MAT-file')
        v = '-v7';
    else
        v = '-v7';
    end
end

function multi_fprintf(fid, varargin)
    fprintf(varargin{:});
    fprintf(fid, varargin{:});
end
