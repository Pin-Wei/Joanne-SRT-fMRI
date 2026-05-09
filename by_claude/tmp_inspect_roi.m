addpath('/home/aclexp/Software/conn');
r = load('/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/no_PM_260421/results/secondlevel/gPPI/str_r12_ran/ROI.mat');
disp('--- top-level fields ---'); disp(fieldnames(r));
if isfield(r,'ROI')
    fprintf('numel ROI: %d\n', numel(r.ROI));
    if numel(r.ROI)>0
        disp('--- ROI(1) fields ---'); disp(fieldnames(r.ROI(1)));
    end
end
% names list
if isfield(r,'names2')
    fprintf('numel names2: %d\n', numel(r.names2));
    disp('--- first 30 names2 ---');
    for i=1:min(30,numel(r.names2)), fprintf('  %d: %s\n', i, r.names2{i}); end
end
if isfield(r,'names')
    fprintf('numel names: %d\n', numel(r.names));
    disp('--- first 30 names ---');
    for i=1:min(30,numel(r.names)), fprintf('  %d: %s\n', i, r.names{i}); end
end
