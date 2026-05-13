function report_conn_out(varargin)
    clc; close all;

    % Configuration -------------------------------------------------------
    
    p = inputParser;
    addOptional(p, 'IP', 23);
    addOptional(p, 'PROJECT', 'no_PM_260512');
    addOptional(p, 'PRESET', 3);
    addOptional(p, 'ROI', 3);
    addOptional(p, 'OVERWRITE', false);
    parse(p, varargin{:});

    IP = p.Results.IP;
    PROJECT = p.Results.PROJECT;
    PRESET = p.Results.PRESET;
    ROI = p.Results.ROI;
    OVERWRITE = p.Results.OVERWRITE;

    if IP == 37, ROOT_DIR = '/media/data3/Joanne_SRT_pw/';
    elseif IP == 23, ROOT_DIR = '/home/aclexp/pinwei/Joanne_SRT_fMRI/'; 
    else, error('Root directory for this IP address has not yet been defined.'); 
    end
    
    PROJ_MAT = fullfile(ROOT_DIR, 'conn_out', [PROJECT, '.mat']);
    PROJ_DIR = fullfile(ROOT_DIR, 'conn_out', PROJECT);

    ANALYSIS = 'gPPI';
    ANA2_DIR = fullfile(PROJ_DIR, 'results', 'secondlevel', ANALYSIS);
    
    ROI_PREFIXS = {'atlas', 'networks', 'joanne'};
    ROI_PREFIX = ROI_PREFIXS{ROI};

    PRESET_NAMES = {'FNC', 'SPC', 'TFCE'};
    PRESET_NAME  = PRESET_NAMES{PRESET};
    
    OUT_DIR   = fullfile(ROOT_DIR, 'reports', PROJECT, ['ROI-', ROI_PREFIX, '_', PRESET_NAME]);
    FIGS_DIR  = fullfile(OUT_DIR, 'figs');
    STATS_DIR = fullfile(OUT_DIR, 'stats');
    BM_DIR    = fullfile(OUT_DIR, 'bookmarks');
    for fd = {FIGS_DIR, STATS_DIR, BM_DIR}, if ~isfolder(fd{1}), mkdir(fd{1}); end; end
    
    assert(exist(PROJ_MAT, 'file'), 'CONN project not exists.');
    load(PROJ_MAT, 'CONN_x');
    CONTRASTS = CONN_x.Results.saved.names;
    
    % Main ----------------------------------------------------------------
    
    summary = cell(numel(CONTRASTS), 7);
    hdr = {'contrast', 'n_clusters', 'ring_jpg', 'matrix_jpg', 'brain_jpg', 'stats_csv', 'note'};
    
    for ic = 1:numel(CONTRASTS)
        contrast = CONTRASTS{ic};
        roiMat = fullfile(ANA2_DIR, contrast, 'ROI.mat');
        n_clusters = NaN; ringFig = ''; matFig = ''; brainFig = ''; statCsv = '';
        
        try
            hfig = conn_displayroi('initfile', roiMat);
            drawnow;        
            data = get(hfig, 'userdata');
    
            % Select ROIs
            roiFilter = make_roi_filter(roiNames, ROI_PREFIX);
            roiNames = data.names(roiFilter);
    
            % Pre-set HC clustering options
            data.clusters_options = struct('type', 'hc', 'groups', NaN, 'param', 0.05);
            set(hfig, 'userdata', data);
    
            conn_displayroi(hfig, [], 'roi.select', [], roiNames);
            drawnow;
    
            % Cluster-level inference
            conn_displayroi(hfig, [], 'fwec.option', PRESET);
            drawnow;
    
            % Fix paper size for consistent JPG aspect (1144 x 624)
            set(hfig, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'PaperPosition', [0 0 11.44 6.24]);
    
            % Save figure 1: ring
            ringFig = fullfile(FIGS_DIR, sprintf('%s (ring).jpg', contrast));
            if ~exist(ringFig, 'file') || OVERWRITE
                conn_displayroi(hfig, [], 'ring_print', ringFig);
                drawnow;
            end
    
            % Save figure 2: matrix 
            matFig = fullfile(FIGS_DIR, sprintf('%s (matrix).jpg', contrast));
            if ~exist(matFig, 'file') || OVERWRITE
                conn_displayroi(hfig, [], 'matrix_print', matFig);
                drawnow;
            end
    
            % Save figure 3: glass brain projections
            brainFig = fullfile(FIGS_DIR, sprintf('%s (brain).jpg', contrast));
            if ~exist(brainFig, 'file') || OVERWRITE
                conn_displayroi(hfig, [], 'glass_print', brainFig);
                drawnow;
            end
    
            % Save bookmarks 
            data = get(hfig, 'userdata');
            bmMat = fullfile(BM_DIR, [contrast, '.mat']);
    
            if ~exist(bmMat, 'file') || OVERWRITE
                connArgs = {'displayroi', 'initfile', data.initfile}; 
                save(bmMat, 'connArgs', '-v7.3');
        
                bmFig = fullfile(BM_DIR, [contrast, '.jpg']);
                conn_print(hfig, bmFig, '-nogui', '-r50', '-nopersistent');
        
                bmTxt = fullfile(BM_DIR, [contrast, '.txt']);
                fid = fopen(bmTxt, 'w');
                desc = {sprintf('%s | roi prefix = %s | HC | preset %s', contrast, ROI_PREFIX, PRESET_NAME)};
                for i = 1:numel(desc), fprintf(fid, '%s\n', desc{i}); end
                fclose(fid);
            end
    
            % Save stats CSV 
            data = get(hfig, 'userdata');
            statCsv = fullfile(STATS_DIR, [contrast, '.csv']);
            n_clusters = write_stats_csv(data, statCsv, ROI_PREFIX);
    
            close all; close all hidden;
            status = 'ok';
            
        catch err
            fprintf('\r\n%s', err.message);
            status = ['ERROR: ', err.message];
        end
    
        summary(ic, :) = {contrast, n_clusters, ringFig, matFig, brainFig, statCsv, status};
    end
    
    % Summary CSV
    summCsv = fullfile(OUT_DIR, 'summary.csv');
    summTbl = array2table(summary, 'VariableNames', hdr);
    writetable(summTbl, summCsv);
end

%% Other functions

function n_clusters = write_stats_csv(data, statCsv, ROI_PREFIX)
% Write structured CSV from data.list2txt.
% Cluster header formats vary by preset:
% - FNC  >> Cluster K/N   F(df1, df2) = val   p_unc  p_FDR
% - SPC  >> Cluster K/N   Mass = val          p_unc  p_FDR  p_FWE
% - TFCE >> Cluster K/N   TFCE = val          p_unc  p_FDR  p_FWE
% Connection rows always: 
% >> Connection roi1 (xyz1)-roi2 (xyz2)  STAT(df) = val  p_unc  p_FDR

    n_clusters = 0;
    statRows = data.list2txt;
    if isempty(statRows), return; end 

    fid = fopen(statCsv, 'w');
    fprintf(fid, 'Cluster_ID,ROI_1,ROI_2,Statistic,p_unc,p_FDR,p_FWE\n');

    for i = 1:numel(statRows)
        row = strtrim(statRows{i});
        if isempty(row), continue; end

        if startsWith(row, 'Cluster')
            roi1 = '--'; roi2 = '--'; 
            parts = regexp(row, '^Cluster\s+(\d+)\/\d+\s+(.*\s*=\s*[\d.]+)\s+(.*)', 'tokens', 'once');
            cid = str2double(parts{1}); 
            statRes = strtrim(parts{2});
            statRes = replace(statRes, ',', ' ');
            pVals = regexp(parts{3}, '\s{2,}', 'split');
            pUnc = strtrim(pVals{1});
            pFDR = strtrim(pVals{2}); 
            if numel(pVals) >= 3, pFWE = strtrim(pVals{3}); else, pFWE = ''; end
            n_clusters = n_clusters + 1;

        else
            parts = regexp(row, '(\\.*\\)\s+(.*\s*=\s*-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)', 'tokens', 'once');
            part1s = split(parts{1}, '\\');
            roi1 = replace(part1s{3}, ',', ' ');
            roi1 = rename_roi(roi1, ROI_PREFIX);
            roi2 = replace(part1s{6}, ',', ' ');
            roi2 = rename_roi(roi2, ROI_PREFIX);
            statRes = strtrim(parts{2});
            pUnc = strtrim(parts{3});
            pFDR = strtrim(parts{4});
            pFWE = '';
        end
    
        fprintf(fid, '%d,%s,%s,%s,%s,%s,%s\n', ...
            cid, roi1, roi2, statRes, pUnc, pFDR, pFWE);
    end

    fclose(fid);
end

function roiFilter = make_roi_filter(roiNames, ROI_PREFIX)
    if matches(ROI_PREFIX, 'joanne')
        roiSetFS = { ...
            'precentral', ...       % SL system
            'Caudate', ...          % SL system
            'Putamen', ...          % SL system
            'Cerebellum', ...       % SL system
            'inferiorparietal', ... % IF system
            'superiorparietal', ... % IF system
            'insula', ...           % IF system
            'CC_Anterior', ...      % IF system
            'Amygdala', ...         % IF system
            'Thalamus' ...          % IF system
        };
        roiSetBA = { ...
            'BA6', ... % supplementary motor area (SMA); SL system
            'BA4', ... % primary motor cortex (M1); SL system
            'BA4' ...  % dorsolateral prefrontal cortex (dlPFC); IF system
        };
        roiFilter = cellfun(@(s) ...
            startsWith(s, 'FS') && any(contains(roiSetFS, s)) || ...
            startsWith(s, 'BA') && any(contains(roiSetBA, s)), ...
            roiNames ...
        );
    else
        roiFilter = cellfun(@(s) startsWith(s, ROI_PREFIX), roiNames);
    end
end

function roi = rename_roi(roi, ROI_PREFIX)
    if matches(ROI_PREFIX, 'joanne')
        roi = replace(roi, ...
            ["ctx-rh-", "R "], ...
            ["ctx-lh-", "L "], ...
            ["BA6", "SMA (Supplementary Motor Area)"], ...
            ["BA4", "M1 (Primary Motor Cortex)"], ...
            ["BA4", "dlPFC (DorsoLateral Prefrontal Cortex)"], ...
            ["precentral", "PreCG (Precentral Gyrus)"], ...
            ["inferiorparietal", "IPL (Inferior Parietal Lobule)"], ...
            ["superiorparietal", "SPL (Superior Parietal Lobule)"], ...
            ["insula", "Insula"], ...
            ["CC_Anterior", "ACC (Anterior Cingulate Cortex)"], ...
            ["-", " "] ...
        );
    else
        roi = replace(roi, [ROI_PREFIX, '.'], ' ');
    end
end