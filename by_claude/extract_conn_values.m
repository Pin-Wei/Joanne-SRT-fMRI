function newNames = extract_conn_values(connFile, analysisName, requests, varargin)
%EXTRACT_CONN_VALUES  Extract per-subject ROI-to-ROI connectivity values from
% CONN 2nd-level analyses and register them as 2nd-level covariates.
%
% Reproduces the CONN GUI "Import values" button in the ROI-level results
% explorer (conn_displayroi.m, callback 'import_values') without launching
% the GUI: reads ROI.mat directly and calls conn_importl2covariate.
%
% USAGE
%   names = extract_conn_values(connFile, analysisName, requests)
%   names = extract_conn_values(..., 'Save', false)
%
% INPUTS
%   connFile     : full path to a CONN project .mat
%   analysisName : 2nd-level analysis name (subfolder under
%                  results/secondlevel/, e.g. 'gPPI', 'SBC_01')
%   requests     : struct array; each element specifies one extraction:
%     .contrast  (required) contrast subfolder under the analysis
%     .source    (required) source ROI name (exact, or unique substring)
%     .target    (required) target ROI name (exact, or unique substring)
%     .condition (optional) condition name (char) or index (numeric); if
%                empty/omitted, every condition in the analysis is exported
%                as a separate covariate.
%     .covName   (optional) base covariate name. With multiple conditions a
%                ' [<condName>]' suffix is appended automatically. Default
%                name: '<contrast> : <source> -> <target> [<cond>]'
%
% OPTIONS (name-value)
%   'Save'    (default true)  call `conn save` after import so covariates
%             persist in the CONN project file.
%   'Verbose' (default true)  print one line per imported covariate.
%
% OUTPUT
%   newNames  cell array of covariate names that were imported.
%
% REQUIREMENTS
%   CONN toolbox on the MATLAB path (uses conn, conn_importl2covariate,
%   conn_loadmatfile). No graphics needed.
%
% NOTE
%   The .mat layout is ROI(i) where i indexes the source ROI; targets are
%   columns of ROI(i).y. Source list = ROI(1).names; target list =
%   ROI(1).names2. Connectivity for subject s, source iS, target iT,
%   condition iC is ROI(iS).y(s, iT, iC).

opts = local_parse_opts(struct('Save',true,'Verbose',true), varargin);

global CONN_x
conn('load', connFile);
analysisRoot = fullfile(CONN_x.folders.secondlevel, analysisName);
if ~isfolder(analysisRoot)
    error('extract_conn_values:NoAnalysis', ...
        'Analysis folder not found: %s', analysisRoot);
end
nsubjects = CONN_x.Setup.nsubjects;

allNames  = {};
allValues = {};

% Cache loaded ROI.mat per contrast to avoid repeated disk I/O.
roiCache = containers.Map;

for k = 1:numel(requests)
    req = requests(k);
    local_require_fields(req, {'contrast','source','target'}, k);

    roiMat = fullfile(analysisRoot, req.contrast, 'ROI.mat');
    if ~isfile(roiMat)
        error('extract_conn_values:NoROImat', ...
            'ROI.mat not found: %s', roiMat);
    end

    if isKey(roiCache, roiMat)
        S = roiCache(roiMat);
    else
        S = conn_loadmatfile(roiMat);   % returns struct with field ROI
        roiCache(roiMat) = S;
    end
    ROI = S.ROI;
    sourceNames = ROI(1).names;     % length = numel(ROI)
    targetNames = ROI(1).names2;    % length = size(ROI(i).y, 2)
    condNames   = ROI(1).ynames;
    selSub      = ROI(1).xX.SelectedSubjects;

    iS = local_lookup(sourceNames, req.source, 'source ROI');
    iT = local_lookup(targetNames, req.target, 'target ROI');

    if isfield(req,'condition') && ~isempty(req.condition)
        iC = local_lookup_cond(condNames, req.condition);
    else
        iC = 1:numel(condNames);
    end

    multiCond = numel(iC) > 1;
    for cc = iC(:)'
        rawY = ROI(iS).y(:, iT, cc);        % [nnz(sel) x 1]
        y = local_expand_subjects(rawY, selSub, nsubjects);

        if isfield(req,'covName') && ~isempty(req.covName)
            nm = req.covName;
            if multiCond, nm = sprintf('%s [%s]', nm, condNames{cc}); end
        else
            nm = sprintf('%s : %s -> %s [%s]', req.contrast, ...
                sourceNames{iS}, targetNames{iT}, condNames{cc});
        end
        allNames{end+1}  = nm;   %#ok<AGROW>
        allValues{end+1} = y;    %#ok<AGROW>
        if opts.Verbose
            fprintf('  + %-90s n=%d (NaN=%d)\n', nm, numel(y), sum(isnan(y)));
        end
    end
end

if isempty(allValues)
    warning('extract_conn_values:Empty', 'Nothing to import.');
    newNames = {};
    return
end

conn_importl2covariate(allNames, allValues, false);
fprintf('Imported %d 2nd-level covariate(s).\n', numel(allNames));

if opts.Save
    conn save
    fprintf('CONN project saved: %s\n', CONN_x.filename);
end

newNames = allNames;
end

% -------------------------------------------------------------------------
function opts = local_parse_opts(opts, args)
for k = 1:2:numel(args)
    fn = args{k};
    if ~isfield(opts, fn)
        error('extract_conn_values:UnknownOpt', 'Unknown option: %s', fn);
    end
    opts.(fn) = args{k+1};
end
end

function local_require_fields(req, fields, k)
for i = 1:numel(fields)
    if ~isfield(req, fields{i}) || isempty(req.(fields{i}))
        error('extract_conn_values:MissingField', ...
            'requests(%d).%s is required', k, fields{i});
    end
end
end

function idx = local_lookup(names, query, label)
% exact first; otherwise unique substring (case-insensitive)
if isnumeric(query), idx = query; return; end
idx = find(strcmp(names, query));
if isempty(idx), idx = find(contains(names, query, 'IgnoreCase', true)); end
if isempty(idx)
    error('extract_conn_values:NotFound', '%s not found: "%s"', label, query);
end
if numel(idx) > 1
    sample = strjoin(names(idx(1:min(end,5))), ', ');
    error('extract_conn_values:Ambiguous', ...
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
        error('extract_conn_values:CondNotFound', ...
            'condition not found: "%s"', query{n});
    end
    if numel(j) > 1
        error('extract_conn_values:CondAmbiguous', ...
            'Ambiguous condition "%s": %d matches', query{n}, numel(j));
    end
    idx(n) = j;
end
end

function y = local_expand_subjects(rawY, selSub, nsubjects)
% rawY : [nnz(selSub) x 1] (single condition slice)
% selSub: logical/0-1 mask of length nsubjects
selSub = logical(selSub(:));
if isempty(selSub) || ~any(selSub)
    y = rawY(:);
    return
end
nnzSel = nnz(selSub);
if rem(numel(rawY), nnzSel) ~= 0
    error('extract_conn_values:SubjectMismatch', ...
        'rawY length %d not a multiple of nnz(SelectedSubjects)=%d', ...
        numel(rawY), nnzSel);
end
rep = numel(rawY) / nnzSel;
y = nan(rep * numel(selSub), 1);
y(repmat(selSub, rep, 1)) = rawY(:);
% Common case (rep==1) yields length == nsubjects.
if rep == 1 && numel(y) ~= nsubjects
    warning('extract_conn_values:SubjectCount', ...
        'expanded length %d != CONN_x.Setup.nsubjects=%d', ...
        numel(y), nsubjects);
end
end
