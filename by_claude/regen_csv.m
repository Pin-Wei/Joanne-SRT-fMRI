% regen_csv.m - regenerate ONLY the structured stats CSV for given contrasts.
% Used to re-parse list2txt with a fixed parser without redoing JPGs/bookmarks.
%
% Reuses write_structured_stats() from run_r2r_networks_bookmarks.m by adding
% that script's directory to the path.

addpath('/home/aclexp/Software/conn');
addpath('/home/aclexp/Software/spm');
addpath('/home/aclexp/pinwei/Joanne_SRT_fMRI/by_claude');

PROJECT  = 'no_PM_260424';
PRESET   = 3;  % TFCE
ANALYSIS = 'gPPI';
CONTRASTS = {'str_r12_ran', 'str_r56_ran'};

PROJECT_ROOT    = '/home/aclexp/pinwei/Joanne_SRT_fMRI';
PROJECT_MAT     = fullfile(PROJECT_ROOT, 'conn_out', [PROJECT '.mat']);
SECONDLEVEL_DIR = fullfile(PROJECT_ROOT, 'conn_out', PROJECT, 'results', 'secondlevel', ANALYSIS);
PRESET_NAMES    = {'FNC','SPC','TFCE'};
preset_tag      = PRESET_NAMES{PRESET};
out_dir         = fullfile(PROJECT_ROOT, 'reports', PROJECT, ['ROI-networks_' preset_tag]);
stats_dir       = fullfile(out_dir, 'stats');
NETWORK_PREFIX  = 'networks.';

fprintf('Loading CONN project...\n');
conn('load', PROJECT_MAT);

for ic = 1:numel(CONTRASTS)
    contrast = CONTRASTS{ic};
    roi_mat  = fullfile(SECONDLEVEL_DIR, contrast, 'ROI.mat');
    fprintf('\n=== %s ===\n', contrast);

    hfig = conn_displayroi('initfile', roi_mat);
    drawnow;
    data = get(hfig, 'userdata');

    is_net = cellfun(@(s) startsWith(s, NETWORK_PREFIX), data.names);
    networks_names = data.names(is_net);
    data.clusters_options = struct('type','hc','groups',NaN,'param',0.05);
    set(hfig, 'userdata', data);

    conn_displayroi(hfig, [], 'roi.select', [], networks_names);
    drawnow;
    conn_displayroi(hfig, [], 'fwec.option', PRESET);
    drawnow;

    data = get(hfig, 'userdata');
    csv_path = fullfile(stats_dir, sprintf('%s.csv', contrast));
    [n_all, n_sig] = write_structured_stats(data, csv_path);
    fprintf('  stats (%d rows, %d sig) -> %s\n', n_all, n_sig, csv_path);

    close(hfig);
end

fprintf('\nDone.\n');
