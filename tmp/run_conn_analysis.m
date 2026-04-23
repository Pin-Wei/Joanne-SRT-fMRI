clear; clc; close all;

%% Setup paths 

IP = 37; 

if IP == 37
    rootDir = '/media/data3/Joanne_SRT_pw/';
    connDir = '/home/aclexp/mytools/matlab/conn';
elseif IP == 23
    rootDir = '/home/aclexp/pinwei/Joanne_SRT_fMRI/';
    connDir = '/home/aclexp/Software/conn';
else 
    error('Root directory for this IP address has not yet been defined.');
end

batchFile = fullfile(rootDir, 'conn_out', 'no_PM_260421.mat');

load(batchFile, 'CONN_x');

%% Configuration

ANALYSIS_NAME = 'gPPI';

COND_OF_INTEREST = [ ...
    "str_r12", "str_r34", "str_r56", "str_r78", ... % "structured"
    "swi_r34", "swi_r56" ... % "switch"
];
COND_LIST = ["random", COND_OF_INTEREST, "incorrect"];

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

BS_COVARS = CONN_x.Setup.l2covariates.names;
if ~contains(BS_COVARS, 'AllSubjects') 
    disp('In order to ');
    exit
end

BS_EFFECT_NAMES = {'AllSubjects'};
BS_CONTRAST = [1];
if contains(BS_COVARS, 'ExcludeOutlierSubjects')
    BS_EFFECT_NAMES = [BS_EFFECT_NAMES, {'ExcludeOutlierSubjects'}];
    BS_CONTRAST = [BS_CONTRAST 0];
end

%% First-level analysis 
% https://web.conn-toolbox.org/fmri-methods/connectivity-measures

if ~contains({CONN_x.Analyses.name}, ANALYSIS_NAME)
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
end

%% Second-level analysis 

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