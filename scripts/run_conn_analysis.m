clear; clc; close all;

%% Setup paths 

IP = 23; % 37

if IP == 37
    rootDir = '/media/data3/Joanne_SRT_pw/';
    connDir = '/home/aclexp/mytools/matlab/conn';
elseif IP == 23
    rootDir = '/home/aclexp/pinwei/Joanne_SRT_fMRI/';
    connDir = '/home/aclexp/Software/conn';
else 
    error('Root directory for this IP address has not yet been defined.');
end

batchFile = fullfile(rootDir, 'conn_out', 'v5_no_param.mat');

%% Configuration

ANALYSIS_NAME = 'gPPI';

COND_OF_INTEREST = [ ...
    "str_r12", "str_r34", "str_r56", "str_r78", ... % "structured"
    "swi_r34", "swi_r56" ... % "switch"
];
COND_LIST = ["random", COND_OF_INTEREST, "incorrect"];

BS_EFFECT_NAMES = {'AllSubjects', 'ExcludeOutlierSubjects'};
BS_CONTRAST = [1 0];

CONTRASTS = struct(); 
for i = 1:numel(COND_OF_INTEREST)
    cond = COND_OF_INTEREST(i);
    CONTRASTS(i).saveas = strcat(cond, '_ran');
    CONTRASTS(i).between_conditions_names = {cond, 'random'};
    CONTRASTS(i).between_conditions_contrast = [1 -1];
end

CONTRASTS(i+1).saveas = 'str_main';
CONTRASTS(i+1).between_conditions_names = {'str_r12', 'str_r34', 'str_r56', 'str_r78', 'random'};
CONTRASTS(i+1).between_conditions_contrast = [1 1 1 1 -4];

CONTRASTS(i+2).saveas = 'swi_main';
CONTRASTS(i+2).between_conditions_names = {'swi_r34', 'swi_r56', 'random'};
CONTRASTS(i+2).between_conditions_contrast = [1 1 -2];

%% First-level analysis 
% https://web.conn-toolbox.org/fmri-methods/connectivity-measures

conn_batch( ...
    'filename', batchFile, ...
    'Analysis.name', ANALYSIS_NAME, ...
    'Analysis.measure', 3, ...     % regression (bivariate)
    'Analysis.modulation', 1, ...  % gPPI interaction effect
    'Analysis.conditions', COND_LIST, ... 
    'Analysis.type', 3, ...        % ROI-to-ROI & Seed-to-Voxel
    'Analysis.done', 1, ...
    'Analysis.overwrite', 'No' ...
);

%% Second-level analysis 

load(batchFile, 'CONN_x');

if ~contains(CONN_x.Setup.l2covariates.names, 'AllSubjects')
    conn_batch( ...
        'filename', batchFile, ...
        'Setup.subjects.groups_names', {'AllSubjects'}, ...
        'Setup.subjects.groups', {ones(CONN_x.Setup.nsubjects, 1)}, ...
        'Setup.subjects.add', 1 ...
    );
end

for i = 1:numel(CONTRASTS)
    try
        C = CONTRASTS(i);
        conn_batch( ...
            'filename', batchFile, ...
            'Results.saveas', C.saveas, ...
            'Results.name', ANALYSIS_NAME, ...
            'Results.between_subjects.effect_names', BS_EFFECT_NAMES, ...
            'Results.between_subjects.contrast', BS_CONTRAST, ...
            'Results.between_conditions.effect_names', C.between_conditions_names, ...
            'Results.between_conditions.contrast', C.between_conditions_contrast, ...
            'Results.display', 0, ...
            'Results.done', 1 ...
        );
    catch err
        fprintf('\nFAILED CONTRAST: %s\n%s\n', C.saveas, err.message);
        continue;
    end
end

%% GLM group-level analyses 
% https://web.conn-toolbox.org/resources/conn-extensions/glm