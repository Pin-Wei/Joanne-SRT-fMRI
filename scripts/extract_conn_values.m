function extract_conn_values(varargin)
    clc; close all;
    
    % Configuration -------------------------------------------------------
    
    p = inputParser;
    addRequired(p, 'IP', @(x) any(x == [23, 37])); % positional arguments
    addRequired(p, 'PROJECT', @ischar);
    addOptional(p, 'ROI', 3, @(x) any(x == [1, 2, 3]));
    parse(p, varargin{:});

    IP = p.Results.IP;
    PROJECT = p.Results.PROJECT;
    ROI = p.Results.ROI;

    ROI_PREFIXS = {'atlas', 'networks', 'joanne'};
    ROI_PREFIX = ROI_PREFIXS{ROI};

    ANALYSIS_NAME = 'gPPI';
    CONTRASTS  = define_for_version(PROJECT);

    % Paths ---------------------------------------------------------------

    if IP == 37, rootDir = '/media/data3/Joanne_SRT_pw/';
    elseif IP == 23, rootDir = '/home/aclexp/pinwei/Joanne_SRT_fMRI/';
    else, error('Root directory for this IP address has not yet been defined.');
    end
    
    resultDir = fullfile(rootDir, 'conn_out', PROJECT, 'results', 'secondlevel', ANALYSIS_NAME);

    outFile = fullfile(rootDir, 'conn_out', 'values', PROJECT, ['ROI_', ROI_PREFIX, '.xlsx']);
    [fd, ~, ~] = fileparts(outFile);
    if ~isfolder(fd), mkdir(fd); end

    % Main ----------------------------------------------------------------

    T = table();

    for ic = 1:numel(CONTRASTS)
        contrast = CONTRASTS{ic};

        roiMat = fullfile(resultDir, contrast, 'ROI.mat');
        if ~isfile(roiMat), error('file not found: %s', roiMat); end

        data = load(roiMat); 
        subjIds = data.summary.design.subjects;
        subjNames = compose('sub_%02d', subjIds);
        nSubjs = numel(subjIds);

        condNames = data.summary.design.conditions;
        nConds = numel(condNames);

        roiAll = data.summary.rois.names;
        roiFilter = make_roi_filter(roiAll,  ROI_PREFIX);
        roiNames = roiAll(roiFilter);
        nRois = numel(roiNames);

        roiData = data.ROI(roiFilter);
        [subjIdxs, roiIdxs, condIdxs] = ndgrid(1:nSubjs, 1:nRois, 1:nConds);

        for iroi = 1:numel(roiData)
            y = roiData(iroi).y(:, roiFilter, :); 
            n = numel(y);

            % fprintf("contrast: %s\n", mat2str(size(repmat({contrast}, n, 1))));
            % fprintf("cond: %s\n", mat2str(size(condNames(condIdxs(:)).')));
            % fprintf("roi1: %s\n", mat2str(size(repmat({roiNames(iroi)}, n, 1))));
            % fprintf("roi2: %s\n", mat2str(size(roiNames(roiIdxs(:)).')));
            % fprintf("subj: %s\n", mat2str(size(subjNames(subjIdxs(:)))));
            % fprintf("y: %s\n", mat2str(size(y(:))));

            tbl = table( ...
                repmat({contrast}, n, 1), ...
                condNames(condIdxs(:)).', ...
                repmat({roiNames(iroi)}, n, 1), ...
                roiNames(roiIdxs(:)).', ...
                subjNames(subjIdxs(:)), ...
                y(:), ...
                'VariableNames', { ...
                    'contrast', ...
                    'condition', ...
                    'roi1', ...
                    'roi2', ...
                    'subject', ...
                    'value' ...
                } ...
            );
            T = [T; tbl];
        end
    end

    T = T(~isnan(T.value), :);
    T.roi1 = cellfun(@(s) rename_roi(s, ROI_PREFIX), T.roi1, 'UniformOutput', false);
    T.roi2 = cellfun(@(s) rename_roi(s, ROI_PREFIX), T.roi2, 'UniformOutput', false);

    writetable(T, outFile);
    fprintf('\r\nFile saved at: %s', outFile);
end

%% Other functions

function CONTRASTS = define_for_version(PROJECT)
    if ismember(PROJECT, {'no_PM_260504', 'no_PM_260515'})
        CONTRASTS = {'str_main', 'swi_main'};
    elseif ismember(PROJECT, {'no_PM_260506', 'no_PM_260512'})
        CONTRASTS = {'str_fst', 'str_snd', 'switch'};
    else
        error("Error! Project not defined.")
        % CONTRASTS = CONN_x.Results.saved.names;
    end
end