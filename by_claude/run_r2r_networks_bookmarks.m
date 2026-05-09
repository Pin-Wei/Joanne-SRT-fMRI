% run_r2r_networks_bookmarks.m
%
% For each ROI-to-ROI second-level contrast in a CONN gPPI project:
%   1. Open second-level display
%   2. Filter to networks.* ROIs
%   3. Apply HC clustering (auto cutoff, lambda=0.05)
%   4. Apply CONN cluster-level inference preset (FNC/SPC/TFCE)
%   5. Save 3 view JPGs: ring / matrix / sagittal
%   6. Save bookmark (txt+mat+jpg)
%   7. Dump structured stats CSV (no \\NNN\\ markers)
%
% Outputs to:
%   reports/<PROJECT>/ROI-networks_<PRESET>/
%       <contrast> (ring).jpg
%       <contrast> (matrix).jpg
%       <contrast> (sagittal view).jpg
%       stats/<contrast>.csv
%       bookmarks/<contrast>_<preset>.{txt,mat,jpg}
%       summary_<PRESET>.csv
%       run_<PRESET>.log
%
% Configure PROJECT and PRESET below; or override via env vars
%   R2R_PROJECT, R2R_PRESET (1=FNC, 2=SPC, 3=TFCE)
%
% Run:
%   matlab -nodesktop -nosplash -batch \
%     "run('/home/aclexp/pinwei/Joanne_SRT_fMRI/by_claude/run_r2r_networks_bookmarks.m')"
%
%   R2R_PROJECT=no_PM_260424 R2R_PRESET=3 matlab -nodesktop -nosplash -batch \
%     "run('/home/aclexp/pinwei/Joanne_SRT_fMRI/by_claude/run_r2r_networks_bookmarks.m')"

addpath('/home/aclexp/Software/conn');
addpath('/home/aclexp/Software/spm');

%% -------- CONFIG ----------------------------------------------------------
PROJECT  = getenv_default('R2R_PROJECT', 'no_PM_260424');
PRESET   = str2double(getenv_default('R2R_PRESET',  '1'));   % 1=FNC, 2=SPC, 3=TFCE
ANALYSIS = 'gPPI';

PROJECT_ROOT    = '/home/aclexp/pinwei/Joanne_SRT_fMRI';
PROJECT_MAT     = fullfile(PROJECT_ROOT, 'conn_out', [PROJECT '.mat']);
PROJECT_FOLDER  = fullfile(PROJECT_ROOT, 'conn_out', PROJECT);
SECONDLEVEL_DIR = fullfile(PROJECT_FOLDER, 'results', 'secondlevel', ANALYSIS);
REPORT_ROOT     = fullfile(PROJECT_ROOT, 'reports', PROJECT);

PRESET_NAMES   = {'FNC', 'SPC', 'TFCE'};
NETWORK_PREFIX = 'networks.';
JPG_DPI        = 300;
%% --------------------------------------------------------------------------

preset_tag   = PRESET_NAMES{PRESET};
out_dir      = fullfile(REPORT_ROOT, ['ROI-networks_' preset_tag]);
fig_dir      = out_dir;
stats_dir    = fullfile(out_dir, 'stats');
bookmark_dir = fullfile(out_dir, 'bookmarks');
log_path     = fullfile(out_dir, 'run.log');

if ~exist(fig_dir,'dir'),      mkdir(fig_dir);      end
if ~exist(stats_dir,'dir'),    mkdir(stats_dir);    end
if ~exist(bookmark_dir,'dir'), mkdir(bookmark_dir); end

fid_log = fopen(log_path, 'w');
logf = @(varargin) multi_fprintf(fid_log, varargin{:});

logf('run_r2r_networks_bookmarks - %s\n', datestr(now));
logf('PROJECT  : %s\n', PROJECT);
logf('PRESET   : %d (%s)\n', PRESET, preset_tag);
logf('ANALYSIS : %s\n', ANALYSIS);
logf('OUT_DIR  : %s\n\n', out_dir);

logf('Loading CONN project...\n');
conn('load', PROJECT_MAT);

d = dir(SECONDLEVEL_DIR);
contrasts = {d([d.isdir] & ~ismember({d.name},{'.','..'})).name};
logf('Found %d contrasts: %s\n\n', numel(contrasts), strjoin(contrasts, ', '));

summary = cell(0,9);
hdr = {'contrast','n_networks','n_clusters','n_sig', ...
       'ring_jpg','matrix_jpg','sagittal_jpg','stats_csv','note'};

for ic = 1:numel(contrasts)
    contrast = contrasts{ic};
    roi_mat  = fullfile(SECONDLEVEL_DIR, contrast, 'ROI.mat');
    logf('\n=== [%d/%d] %s ===\n', ic, numel(contrasts), contrast);

    if ~isfile(roi_mat)
        logf('  [SKIP] ROI.mat missing\n');
        summary(end+1,:) = {contrast, NaN, NaN, NaN, '','','','', 'ROI.mat missing'}; %#ok<SAGROW>
        continue;
    end

    hfig = [];
    try
        hfig = conn_displayroi('initfile', roi_mat);
        drawnow;
        data = get(hfig, 'userdata');

        % --- Filter to networks.* ---
        is_net = cellfun(@(s) startsWith(s, NETWORK_PREFIX), data.names);
        networks_names = data.names(is_net);
        logf('  networks ROIs: %d (of %d total)\n', numel(networks_names), numel(data.names));
        if numel(networks_names) < 2
            error('Fewer than 2 networks ROIs found.');
        end
        % Pre-set HC clustering opts so roi.select doesn't pop interactive dialog
        data.clusters_options = struct('type','hc','groups',NaN,'param',0.05);
        set(hfig, 'userdata', data);

        conn_displayroi(hfig, [], 'roi.select', [], networks_names);
        drawnow;

        % --- Apply preset (cluster-level inference) ---
        conn_displayroi(hfig, [], 'fwec.option', PRESET);
        drawnow;

        % Fix paper size for consistent JPG aspect (matches GUI default 1144x624)
        set(hfig, 'PaperPositionMode','manual', 'PaperUnits','inches', ...
                  'PaperPosition',[0 0 11.44 6.24]);

        % --- View 1: ring (default, no brain) ---
        ring_jpg = fullfile(fig_dir, sprintf('%s (ring).jpg', contrast));
        conn_displayroi(hfig, [], 'view-ring');
        drawnow;
        conn_print(hfig, ring_jpg, '-nogui', sprintf('-r%d', JPG_DPI), '-nopersistent');
        logf('  ring JPG     -> %s\n', ring_jpg);

        % --- View 2: matrix display (separate figure via conn_montage_display) ---
        matrix_jpg = fullfile(fig_dir, sprintf('%s (matrix).jpg', contrast));
        conn_displayroi(hfig, [], 'matrix_print', matrix_jpg);
        drawnow;
        logf('  matrix JPG   -> %s\n', matrix_jpg);

        % --- View 3: glass brain projection (multi-view, no_PM_260421-style)
        % glass_print routes through conn_mesh_display, which spawns its
        % own figure and writes the JPG via varargin{1}. No conn_print needed.
        sagittal_jpg = fullfile(fig_dir, sprintf('%s (sagittal view).jpg', contrast));
        conn_displayroi(hfig, [], 'glass_print', sagittal_jpg);
        drawnow;
        logf('  glass JPG    -> %s\n', sagittal_jpg);
        % Close any leftover hidden mesh-display figures
        all_figs = findall(0, 'type','figure');
        for fi = 1:numel(all_figs)
            if all_figs(fi) ~= hfig, close(all_figs(fi)); end
        end

        % --- Bookmark ---
        stamp   = sprintf('%s_networks_%s_%s', datestr(now,'yyyy_mm_dd_HHMMSSFFF'), preset_tag, contrast);
        bk_base = fullfile(bookmark_dir, [stamp '.bookmark']);
        data    = get(hfig, 'userdata');
        conn_args = {'displayroi','initfile',data.initfile};                 %#ok<NASGU>
        descr   = {sprintf('%s | networks-only | HC | preset %s', contrast, preset_tag)};
        ftxt = fopen([bk_base '.txt'],'w');
        for i = 1:numel(descr), fprintf(ftxt, '%s\n', descr{i}); end
        fclose(ftxt);
        save([bk_base '.mat'], 'conn_args', '-v7.3');
        conn_print(hfig, [bk_base '.jpg'], '-nogui', '-r50', '-nopersistent');
        logf('  bookmark     -> %s.{txt,mat,jpg}\n', bk_base);

        % --- Structured stats CSV ---
        data = get(hfig, 'userdata');
        csv_path = fullfile(stats_dir, sprintf('%s.csv', contrast));
        [n_all, n_sig] = write_structured_stats(data, csv_path);
        logf('  stats (%d rows, %d sig) -> %s\n', n_all, n_sig, csv_path);

        n_clusters = 0;
        if isfield(data,'CLUSTER_labels') && ~isempty(data.CLUSTER_labels)
            n_clusters = max(data.CLUSTER_labels(:));
        elseif isfield(data,'clusters') && ~isempty(data.clusters)
            n_clusters = max(data.clusters);
        end

        summary(end+1,:) = {contrast, numel(networks_names), n_clusters, n_sig, ...
                            ring_jpg, matrix_jpg, sagittal_jpg, csv_path, 'ok'}; %#ok<SAGROW>
        close(hfig); hfig = [];
        % conn_montage_display may leak hidden figures
        close all hidden;
    catch ME
        logf('  [ERROR] %s\n', ME.message);
        for s = 1:min(numel(ME.stack),3)
            logf('           at %s:%d\n', ME.stack(s).name, ME.stack(s).line);
        end
        if ~isempty(hfig) && ishandle(hfig), close(hfig); end
        summary(end+1,:) = {contrast, NaN, NaN, NaN, '','','','', ['ERROR: ' ME.message]}; %#ok<SAGROW>
    end
end

% Summary log
logf('\n\n===== SUMMARY =====\n');
logf('%-16s %5s %5s %5s   %s\n', 'contrast','#net','#clu','#sig','note');
for i = 1:size(summary,1)
    logf('%-16s %5s %5s %5s   %s\n', summary{i,1}, num2str(summary{i,2}), ...
         num2str(summary{i,3}), num2str(summary{i,4}), summary{i,9});
end

% Summary CSV
sfcsv = fopen(fullfile(out_dir, 'summary.csv'), 'w');
fprintf(sfcsv, '%s\n', strjoin(hdr, ','));
for i = 1:size(summary,1)
    fprintf(sfcsv, '%s,%s,%s,%s,%s,%s,%s,%s,"%s"\n', ...
        summary{i,1}, num2str(summary{i,2}), num2str(summary{i,3}), num2str(summary{i,4}), ...
        summary{i,5}, summary{i,6}, summary{i,7}, summary{i,8}, summary{i,9});
end
fclose(sfcsv);

fclose(fid_log);
fprintf('\nDone. Log: %s\n', log_path);
fprintf('Output: %s\n', out_dir);


% --- Helpers ----------------------------------------------------------------

function v = getenv_default(name, default_val)
    v = getenv(name);
    if isempty(v), v = default_val; end
end

function multi_fprintf(fid, varargin)
    fprintf(varargin{:});
    if fid > 0, fprintf(fid, varargin{:}); end
end

% write_structured_stats() now lives in its own file (write_structured_stats.m)
% so that auxiliary scripts (e.g. regen_csv.m) can call it via addpath.
