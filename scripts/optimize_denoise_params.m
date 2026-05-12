function optimize_denoise_params(varargin)
    clc; close all;
    
    % Configuration -------------------------------------------------------
    
    p = inputParser;
    addOptional(p, 'IP', 37);
    addOptional(p, 'PROJECT', 'no_PM_260421');
    addOptional(p, 'PREV_XLS', {}); % paths to existing results (.xlsx)
    parse(p, varargin{:});

    IP = p.Results.IP;
    PROJECT = p.Results.PROJECT;
    PREV_XLS = p.Results.PREV_XLS;

    POLY_ORD_OPTS  = [1, 2, 3];                     % options of polynomial detrending order
    BP_LO          = 0.008;
    BP_HI_OPTS     = [0.09, 999];                   % band-pass / high-pass filter 
    SIMULT_OPTS    = [1, 0];                        % whether to use simultaneous regression & band-pass
    MOT24_OPTS     = [1, 0];                        % whether to add quadratic motion parameters 
    N_ACOMP_OPTS   = [5, 10, 15, 20, 30, 40, 50];   % numbers of aCompCor components
    N_ACOMP_COARSE = [5, 20, 50];
    N_AROMA_OPTS   = [0, 5, 10, 999];               % numbers of ICA-AROMA components
    N_AROMA_COARSE = [0, 999];
    RM_GMR_OPTS    = [1, 0];                        % whether to remove the average signal within the grey matter (GM) mask
    
    MAX_ROWS = 5000;

    % Paths ---------------------------------------------------------------
    
    if IP == 37, rootDir = '/media/data3/Joanne_SRT_pw/';
    elseif IP == 23, rootDir = '/home/aclexp/pinwei/Joanne_SRT_fMRI/';
    else, error('Root directory for this IP address has not yet been defined.');
    end
    
    batchFile = fullfile(rootDir, 'conn_out', [PROJECT, '.mat']);
    assert(isfile(batchFile), sprintf('\r\n%s not found.', batchFile));
    
    t = datetime('now', 'Format', 'yyyy-MM-dd');
    xlsFile = fullfile(rootDir, 'conn_out', sprintf('opt_results_%s.xlsx', t));
    logFile = fullfile(rootDir, 'logs', sprintf('opt_log_%s.log', t));
    
    [fd, ~, ~] = fileparts(logFile);
    if ~isfolder(fd), mkdir(fd); end    

    % Initaliztion --------------------------------------------------------
    
    load(batchFile, 'CONN_x');
    origPreproc = CONN_x.Preproc;
    
    targets = {'QC_InvalidScans', 'QC_ProportionValidScans', 'QC_MeanMotion'};
    qcVars = intersect(targets, CONN_x.Setup.l2covariates.names, 'stable');
    
    paramColumns = {'POLY_ORD', 'BP_LO', 'BP_HI', 'SIMULT', 'MOT24', 'N_ACOMP','N_AROMA', 'RM_GMR'};
    scoreColumns = {'Validity', 'Quality', 'Sensitivity', 'Mean_QC'};
    requiredCols = [paramColumns, scoreColumns];
    resultColumns = [{'Iter'}, {'Stage'}, requiredCols, {'Runtime_s'},];
    R = nan(MAX_ROWS, numel(resultColumns)); % pre-allocate memory
    S = cell(MAX_ROWS, 1);
    
    scoreIdx = find(strcmp(resultColumns, 'Mean_QC'), 1);
    polyIdx  = find(strcmp(resultColumns, 'POLY_ORD'), 1);
    bphIdx   = find(strcmp(resultColumns, 'BP_HI'), 1);
    simIdx   = find(strcmp(resultColumns, 'SIMULT'), 1);
    m24Idx   = find(strcmp(resultColumns, 'MOT24'), 1);
    acompIdx = find(strcmp(resultColumns, 'N_ACOMP'), 1);
    aromaIdx = find(strcmp(resultColumns, 'N_AROMA'), 1);
    gmrIdx   = find(strcmp(resultColumns, 'RM_GMR'), 1);
    
    iter = 0; % running counter
    
    % Start logging -------------------------------------------------------
    
    logFID = fopen(logFile, 'a'); 
    
    t = datetime('now', 'Format', 'MMMM d, yyyy h:mm a'); 
    fprintf(logFID, '\r\n======================== %s ========================', t);
    
    [~, fn, ~] = fileparts(batchFile);
    fprintf(logFID, '\r\nRun optimization on "%s"', fn);
    fprintf(logFID, '\r\n');

    [prevKeys, prevT] = load_prior_results(PREV_XLS, requiredCols, logFID);
    write_log(logFID, '%d tested parameter combinations are loaded.', height(prevT));

    % Stage 1: coarse grid ------------------------------------------------
    
    [g1, g2, g3, g4, g5, g6, g7] = ndgrid( ...
        1:numel(POLY_ORD_OPTS), ...
        1:numel(BP_HI_OPTS), ...
        1:numel(SIMULT_OPTS), ...
        1:numel(MOT24_OPTS), ... 
        1:numel(N_ACOMP_COARSE), ...
        1:numel(N_AROMA_COARSE), ...
        1:numel(RM_GMR_OPTS) ...
    );
    grid = [g1(:), g2(:), g3(:), g4(:), g5(:), g6(:), g7(:)];

    write_log(logFID, '\r\n');
    write_log(logFID, 'Stage 1 starts');
    
    t0 = tic;
    for ic = 1:size(grid, 1)
        col1 = NaN; elapsed = NaN;
        
        poly  = POLY_ORD_OPTS(grid(ic, 1));
        bpHi  = BP_HI_OPTS(grid(ic, 2));
        sim   = SIMULT_OPTS(grid(ic, 3));
        m24   = MOT24_OPTS(grid(ic, 4));
        acomp = N_ACOMP_COARSE(grid(ic, 5));
        aroma = N_AROMA_COARSE(grid(ic, 6));
        gmr   = RM_GMR_OPTS(grid(ic, 7));
    
        key = make_param_key( ...
            poly, BP_LO, bpHi, sim, m24, acomp, aroma, gmr);

        if isKey(prevKeys, key)
            [dv, dq, ds, mean_qc] = find_prior_results( ...
                prevT, poly, BP_LO, bpHi, sim, m24, acomp, aroma, gmr);

            status = 'prev';
        else
            if bpHi == 999, bp = [BP_LO, Inf]; 
            else, bp = [BP_LO, bpHi]; 
            end

            [dv, dq, ds, mean_qc, elapsed, status] = run_one_combi( ...
                origPreproc, poly, bp, sim, m24, acomp, aroma, gmr, qcVars);

            prevKeys(key) = true; 
            iter = iter + 1;
            col1 = iter;
        end

        R(ic, :) = [col1, 1, ...
            poly, BP_LO, bpHi, sim, m24, acomp, aroma, gmr, ...
            dv, dq, ds, mean_qc, elapsed];

        S{ic} = status;
    end
    
    write_log(logFID, 'Stage 1 completed (%d iters, Elapsed: %.1f min)', iter, toc(t0)/60);  
    nIterS1 = iter;
    offset = ic;
    
    % Find best iter ------------------------------------------------------
    
    [~, bestIdx] = max(R(:, scoreIdx));
    bestPoly  = R(bestIdx, polyIdx);
    bestBp2   = R(bestIdx, bphIdx);
    bestSim   = R(bestIdx, simIdx);
    bestM24   = R(bestIdx, m24Idx);
    bestAcomp = R(bestIdx, acompIdx);
    bestAroma = R(bestIdx, aromaIdx);
    bestGmr   = R(bestIdx, gmrIdx);
    
    write_log(logFID, 'Best iter (%d): POLY=%d  BP=[%.3f, %.3g]  SIM=%s  M24=%s  ACOMP=%d  AROMA=%g  GMR=%s', ...
        bestIdx, bestPoly, BP_LO, bestBp2, log2str(bestSim), log2str(bestM24), bestAcomp, bestAroma, log2str(bestGmr));
    
    % Stage 2: fine grid around winner ------------------------------------
    
    [g1, g2] = ndgrid( ...
        1:numel(N_ACOMP_OPTS), ...
        1:numel(N_AROMA_OPTS) ...
    );
    grid2 = [g1(:), g2(:)];
    
    write_log(logFID, '\r\n');
    write_log(logFID, 'Stage 2 starts');
    
    t0 = tic;
    for ic = 1:size(grid2, 1)
        rowIdx = offset + ic;
        dv = NaN; dq = NaN; ds = NaN; mean_qc = NaN; col1 = NaN; elapsed = NaN;

        acomp = N_ACOMP_OPTS(grid2(ic, 1));
        aroma = N_AROMA_OPTS(grid2(ic, 2));
    
        key = make_param_key( ...
            bestPoly, BP_LO, bestBp2, bestSim, bestM24, acomp, aroma, bestGmr);

        if isKey(prevKeys, key)
            try
                [dv, dq, ds, mean_qc] = find_prior_results( ...
                    prevT, bestPoly, BP_LO, estBp2, bestSim, bestM24, acomp, aroma, bestGmr);

                status = 'Prev';

            catch err
                status = ['FAIL: ', strrep(err.message, newline, ' ')];
            end

        else
            if bestBp2 == 999, bestBp = [BP_LO, Inf]; 
            else, bestBp = [BP_LO, bestBp2]; 
            end

            [dv, dq, ds, mean_qc, elapsed, status] = run_one_combi( ...
                origPreproc, bestPoly, bestBp, bestSim, bestM24, acomp, aroma, bestGmr, qcVars);
        
            prevKeys(key) = true; 
            iter = iter + 1;
            col1 = iter;
        end

        R(rowIdx, :) = [col1, 2, ...
            bestPoly, BP_LO, bestBp2, bestSim, bestM24, acomp, aroma, bestGmr, ...
            dv, dq, ds, mean_qc, elapsed];

        S{rowIdx} = status;
    end
    
    write_log(logFID, 'Stage 2 completed (%d iters, Elapsed: %.1f min)', iter - nIterS1, toc(t0)/60);  
    write_log(logFID, 'All %d iterations completed.', iter);
    
    % Write results to file -----------------------------------------------
    
    T = array2table(R(1:rowIdx, :), 'VariableNames', resultColumns);
    T.Status = S(1:rowIdx);
    
    [~, sort_idx] = sort(T.Mean_QC, 'descend', 'MissingPlacement', 'last');
    T = T(sort_idx, :);
    
    writetable(T, xlsFile);
end

%% Other functions

function [dv, dq, ds, mean_qc] = find_prior_results(prevT, poly, bp1, bp2, sim, m24, acomp, aroma, gmr)
    idx = find( ...
        (prevT.POLY_ORD == poly) & ...
        (prevT.BP_LO == bp1) & ...
        (prevT.BP_HI == bp2) & ...
        (prevT.SIMULT == sim) & ...
        (prevT.MOT24 == m24) & ...
        (prevT.N_ACOMP == acomp) & ...
        (prevT.N_AROMA == aroma) & ...
        (prevT.RM_GMR == gmr), ...
        1 ...
    );
    row = prevT(idx, :);
    dv = row.Validity;
    dq = row.Quality;
    ds = row.Sensitivity;
    mean_qc = row.Mean_QC;
end

function [prevKeys, prevT] = load_prior_results(PREV_XLS, requiredCols, logFID)
    prevT = table();
    prevKeys = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    if isempty(PREV_XLS)
        write_log(logFID, 'No existing results loaded.');
        return;
    end

    for ifile = 1:numel(PREV_XLS)
        fp = PREV_XLS{ifile};
        assert(isfile(fp), sprintf('ERROR! File "%s" not exists.', fp));

        try
            T = readtable(fp); 
        catch err
            write_log(logFID, 'WARNNING! Failed to load "%s": %s', fp, err.message);
            continue; 
        end

        assert(all(ismember(requiredCols, T.Properties.VariableNames)), 'Columns not match.');
        write_log(logFID, 'Loading %d results from "%s"', height(T), fp);

        if(isempty(prevT))
            nPrev = 0;
            prevT = T(:, requiredCols);
        else
            nPrev = height(prevT);
            prevT = [prevT(:, requiredCols); T(:, requiredCols)];
            prevT = unique(prevT(:, requiredCols), 'stable');
        end
        
        for irow = 1:height(T)
            k = make_param_key( ...
                T.POLY_ORD(irow), ...
                T.BP_LO(irow), ...
                T.BP_HI(irow), ... 
                T.SIMULT(irow), ... 
                T.MOT24(irow), ... 
                T.N_ACOMP(irow), ... 
                T.N_AROMA(irow), ... 
                T.RM_GMR(irow) ...
            );
            prevKeys(k) = true;
        end
    end
end

function k = make_param_key(poly, bp1, bp2, sim, m24, acomp, aroma, gmr)
    k = sprintf('%d|%d|%d|%s|%s|%d|%d|%s', ...
        poly, ...
        bp1, ...
        bp2, ...
        log2str(sim), ...
        log2str(m24), ...
        acomp, ...
        aroma, ...
        log2str(gmr) ...
    );
end

function [dv, dq, ds, mean_qc, elapsed, status] = run_one_combi( ...
    origPreproc, poly, bp, sim, m24, acomp, aroma, gmr, qcVars ...
)
    global CONN_x;
    dv = NaN; dq = NaN; ds = NaN; mean_qc = NaN; 
    try
        t0 = tic;
        CONN_x.Preproc = origPreproc;
        CONN_x.Preproc.detrending = poly;
        CONN_x.Preproc.filter = bp;
        CONN_x.Preproc.regbp = 1 + sim;
        CONN_x.Preproc.confounds = update_confounds( ...
            origPreproc.confounds, m24, acomp, aroma, gmr);
        dv = safe_score(@() conn_qascores('DataValidity', [], []));
        dq = safe_score(@() conn_qascores('DataQuality', [], [], qcVars, {}));
        ds = safe_score(@() conn_qascores('DataSensitivity', [], [], [], [], 'extreme'));
        mean_qc = mean([dv, dq, ds], 'omitnan');
        elapsed = toc(t0);
        status = 'OK';

    catch err
        elapsed = toc(t0);
        status = ['FAIL: ', strrep(err.message, newline, ' ')];
    end
end

function out = safe_score(fn)
    out = fn(); if isempty(out), out = NaN; end
end

function conf = update_confounds(conf_orig, m24, acomp, aroma, gmr)
    conf = conf_orig;

    idx = find(strcmp(conf.names, 'realignment'), 1);
    assert(~isempty(idx), '"realignment" not in confounds.names');
    conf.deriv{idx} = 1; % always 1st derivative
    if m24, conf.power{idx} = 2; else, conf.power{idx} = 1; end 

    idx = find(strcmp(conf.names, 'aCompCor'), 1);
    assert(~isempty(idx), '"aCompCor" not in confounds.names');
    conf.dimensions{idx} = acomp;

    idx = find(strcmp(conf.names, 'aroma'), 1);
    if aroma == 0 && ~isempty(idx)
        conf = remove_confound(conf, idx);
    elseif isempty(idx)
        conf = append_confound(conf, 'aroma', 'cov', aroma);
    else
        conf.dimensions{idx} = aroma;
    end

    idx = find(strcmp(conf.names, 'Grey Matter'), 1);
    if gmr && isempty(idx)
        conf = append_confound(conf, 'Grey Matter', 'roi', 1);
    elseif ~gmr && ~isempty(idx)
        conf = remove_confound(conf, idx);
    end
end

function conf = append_confound(conf, name, type, dim)
    conf.names{end+1} = name;
    conf.types{end+1} = type;
    conf.power{end+1} = 1;
    conf.deriv{end+1} = 0;
    conf.dimensions{end+1} = [dim, dim];
    if isfield(conf,'filter'), conf.filter{end+1} = 0; end
    if isfield(conf,'fixed'), conf.fixed{end+1} = 0; end
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

function write_log(fid, fmt, varargin)
    msg  = sprintf(fmt, varargin{:});
    t = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:SS');
    line = sprintf('\r\n[%s] %s', t, msg);
    fprintf('%s', line);
    if fid > 0, fprintf(fid, '%s', line); end
end

function str = log2str(a)
    if a, str = 'true'; else, str = 'false'; end
end