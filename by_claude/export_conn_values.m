function [T, csvPath] = export_conn_values(connFile, analysisName, requests, varargin)
%EXPORT_CONN_VALUES  Export per-subject ROI-to-ROI connectivity values to a
% long-format CSV table (does NOT register them as CONN 2nd-level covariates;
% use extract_conn_values.m for that).
%
% Reads ROI.mat directly from a CONN 2nd-level analysis folder; no GUI needed.
%
% USAGE
%   [T, csv] = export_conn_values(connFile, analysisName, requests)
%   [T, csv] = export_conn_values(..., 'OutFile', '/path/out.csv')
%   [T, csv] = export_conn_values(..., 'DropNaN', true)
%
% INPUTS
%   connFile     full path to a CONN project .mat
%   analysisName 2nd-level analysis name (subfolder under
%                results/secondlevel/, e.g. 'gPPI')
%   requests     struct array; each element specifies one extraction:
%     .contrast  (required) contrast subfolder under the analysis
%     .source    (required) source ROI name (exact, or unique substring)
%     .target    (required) target ROI name (exact, or unique substring)
%     .condition (optional) condition name or index. If empty/omitted, every
%                condition is exported as its own rows.
%
% OPTIONS (name-value)
%   'OutFile'  default: by_claude/output/conn_values_<analysis>_<stamp>.csv
%             (set [] or '' to skip writing the CSV — only returns T)
%   'DropNaN'  default false. If true, omit rows for subjects with NaN value
%             (i.e. subjects deselected via SelectedSubjects).
%   'Verbose'  default true.
%
% OUTPUT
%   T        MATLAB table in LONG format with columns:
%            subject, contrast, source, target, condition, value
%   csvPath  full path of the written CSV (empty if 'OutFile' was [] or '').

opts = local_parse_opts(struct( ...
    'OutFile', '__auto__', ...
    'DropNaN', false, ...
    'Verbose', true), varargin);

global CONN_x
conn('load', connFile);
analysisRoot = fullfile(CONN_x.folders.secondlevel, analysisName);
if ~isfolder(analysisRoot)
    error('export_conn_values:NoAnalysis', ...
        'Analysis folder not found: %s', analysisRoot);
end
nsubjects = CONN_x.Setup.nsubjects;

subject   = [];
contrast  = {};
source    = {};
target    = {};
condition = {};
value     = [];

roiCache = containers.Map;

for k = 1:numel(requests)
    req = requests(k);
    local_require_fields(req, {'contrast','source','target'}, k);

    roiMat = fullfile(analysisRoot, req.contrast, 'ROI.mat');
    if ~isfile(roiMat)
        error('export_conn_values:NoROImat', 'ROI.mat not found: %s', roiMat);
    end

    if isKey(roiCache, roiMat)
        S = roiCache(roiMat);
    else
        S = conn_loadmatfile(roiMat);
        roiCache(roiMat) = S;
    end
    ROI = S.ROI;
    sourceNames = ROI(1).names;
    targetNames = ROI(1).names2;
    condNames   = ROI(1).ynames;
    selSub      = ROI(1).xX.SelectedSubjects;

    iS = local_lookup(sourceNames, req.source, 'source ROI');
    iT = local_lookup(targetNames, req.target, 'target ROI');

    if isfield(req,'condition') && ~isempty(req.condition)
        iC = local_lookup_cond(condNames, req.condition);
    else
        iC = 1:numel(condNames);
    end

    for cc = iC(:)'
        rawY = ROI(iS).y(:, iT, cc);
        y = local_expand_subjects(rawY, selSub, nsubjects);
        n = numel(y);

        subject   = [subject;   (1:n).'];                        %#ok<AGROW>
        contrast  = [contrast;  repmat({req.contrast},       n,1)]; %#ok<AGROW>
        source    = [source;    repmat({sourceNames{iS}},    n,1)]; %#ok<AGROW>
        target    = [target;    repmat({targetNames{iT}},    n,1)]; %#ok<AGROW>
        condition = [condition; repmat({condNames{cc}},      n,1)]; %#ok<AGROW>
        value     = [value;     y];                              %#ok<AGROW>

        if opts.Verbose
            fprintf('  + %s : %s -> %s [%s]   n=%d (NaN=%d)\n', ...
                req.contrast, sourceNames{iS}, targetNames{iT}, ...
                condNames{cc}, n, sum(isnan(y)));
        end
    end
end

T = table(subject, contrast, source, target, condition, value);

if opts.DropNaN
    T = T(~isnan(T.value), :);
end

csvPath = '';
if ~(isempty(opts.OutFile) && ~ischar(opts.OutFile))  % allow [] to skip
    if ischar(opts.OutFile) && strcmp(opts.OutFile, '__auto__')
        outDir = fullfile(fileparts(mfilename('fullpath')), 'output');
        if ~isfolder(outDir), mkdir(outDir); end
        stamp = datestr(now, 'yyyymmdd_HHMMSS');
        csvPath = fullfile(outDir, ...
            sprintf('conn_values_%s_%s.csv', analysisName, stamp));
    elseif ischar(opts.OutFile) && ~isempty(opts.OutFile)
        csvPath = opts.OutFile;
        outDir = fileparts(csvPath);
        if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
    else
        csvPath = '';
    end
    if ~isempty(csvPath)
        writetable(T, csvPath);
        if opts.Verbose
            fprintf('Wrote %d rows -> %s\n', height(T), csvPath);
        end
    end
end
end

% -------------------------------------------------------------------------
function opts = local_parse_opts(opts, args)
for k = 1:2:numel(args)
    fn = args{k};
    if ~isfield(opts, fn)
        error('export_conn_values:UnknownOpt', 'Unknown option: %s', fn);
    end
    opts.(fn) = args{k+1};
end
end

function local_require_fields(req, fields, k)
for i = 1:numel(fields)
    if ~isfield(req, fields{i}) || isempty(req.(fields{i}))
        error('export_conn_values:MissingField', ...
            'requests(%d).%s is required', k, fields{i});
    end
end
end

function idx = local_lookup(names, query, label)
if isnumeric(query), idx = query; return; end
idx = find(strcmp(names, query));
if isempty(idx), idx = find(contains(names, query, 'IgnoreCase', true)); end
if isempty(idx)
    error('export_conn_values:NotFound', '%s not found: "%s"', label, query);
end
if numel(idx) > 1
    sample = strjoin(names(idx(1:min(end,5))), ', ');
    error('export_conn_values:Ambiguous', ...
        'Ambiguous %s "%s": %d matches (e.g. %s)', label, query, numel(idx), sample);
end
end

function idx = local_lookup_cond(condNames, query)
if isnumeric(query), idx = query; return; end
if ischar(query), query = {query}; end
idx = zeros(1, numel(query));
for n = 1:numel(query)
    j = find(strcmp(condNames, query{n}));
    if isempty(j), j = find(contains(condNames, query{n}, 'IgnoreCase', true)); end
    if isempty(j)
        error('export_conn_values:CondNotFound', ...
            'condition not found: "%s"', query{n});
    end
    if numel(j) > 1
        error('export_conn_values:CondAmbiguous', ...
            'Ambiguous condition "%s": %d matches', query{n}, numel(j));
    end
    idx(n) = j;
end
end

function y = local_expand_subjects(rawY, selSub, nsubjects)
selSub = logical(selSub(:));
if isempty(selSub) || ~any(selSub)
    y = rawY(:);
    return
end
nnzSel = nnz(selSub);
if rem(numel(rawY), nnzSel) ~= 0
    error('export_conn_values:SubjectMismatch', ...
        'rawY length %d not a multiple of nnz(SelectedSubjects)=%d', ...
        numel(rawY), nnzSel);
end
rep = numel(rawY) / nnzSel;
y = nan(rep * numel(selSub), 1);
y(repmat(selSub, rep, 1)) = rawY(:);
if rep == 1 && numel(y) ~= nsubjects
    warning('export_conn_values:SubjectCount', ...
        'expanded length %d != CONN_x.Setup.nsubjects=%d', numel(y), nsubjects);
end
end
