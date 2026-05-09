% test_conn_save.m
%
% 目的：驗證在兩個位置補上 `conn save` 之後，
%   (1) ExcludeOutlierSubjects 真的會寫進 CONN_x.Setup.l2covariates
%   (2) Results.saveas 真的會寫進 CONN_x.Results.saved
%
% 規模刻意縮到最小，但保留 ROI-to-ROI 分析所需的多 seed 結構：
%   2 名受試者 × 2 個 run × 3 個 condition × 3 個 seed (DefaultMode 網路) × 1 個 contrast
%
% Setup 載入完整的 networks.nii（multi-label，32 ROIs），但 Analysis.sources
% 只指定 3 顆 DefaultMode 節點作為 seed，第一層分析就只跑這 3 顆。
%
% 跑法：
%   cd /home/aclexp/pinwei/Joanne_SRT_fMRI/by_claude
%   matlab -batch "addpath('/home/aclexp/Software/conn'); test_conn_save"

clear; clc; close all;

%% 1. Configuration

PROJECT  = 'tmp';
rootDir  = '/home/aclexp/pinwei/Joanne_SRT_fMRI/';
connDir  = '/home/aclexp/Software/conn';
bidsDir  = fullfile(rootDir, 'data', 'fmriprep');

batchFile = fullfile(rootDir, 'conn_out', [PROJECT, '.mat']);
logFile   = fullfile(rootDir, 'logs',     [PROJECT, '.log']);

[fd, ~, ~] = fileparts(logFile);   if ~isfolder(fd), mkdir(fd); end
[fd, ~, ~] = fileparts(batchFile); if ~isfolder(fd), mkdir(fd); end

% 每次跑都從零開始
projDir = fullfile(rootDir, 'conn_out', PROJECT);
if exist(batchFile, 'file'), delete(batchFile); end
if isfolder(projDir),        rmdir(projDir, 's'); end

SID_LIST = [1, 2];
N_SUBJS  = numel(SID_LIST);
N_RUNS   = 2;
TR       = 2.0;
SPACE    = 'MNI152NLin2009cAsym';

% gPPI 至少要兩個 condition，這裡用三個
COND_NAMES = ["random", "str_fst_r1", "str_fst_r2"];
N_CONDS    = numel(COND_NAMES);

COVAR_NAMES = {'realignment', 'fd', 'scrubbing', 'aCompCor', 'aroma'};
N_COVARS    = numel(COVAR_NAMES);

L2_COVARS = {'QC_MeanMotion', 'QC_InvalidScans'};

% Denoising
FWHM     = 6;
POLY_ORD = 3;
BP_HZ    = [0.008 0.09];
N_ACOMP  = 5;
N_AROMA  = 10;
RM_GMR   = true;
MOT24    = true;

CONFOUND_NAMES = {'realignment', 'scrubbing', 'aCompCor', 'aroma'};
CONFOUND_NAMES = [CONFOUND_NAMES, append('Effect of ', COND_NAMES)];
if RM_GMR, CONFOUND_NAMES = ['Grey Matter', CONFOUND_NAMES]; end

% Setup 載入完整 networks 圖譜（32 ROIs）
ROI_NAMES = {'networks'};
ROI_FILES = {fullfile(connDir, 'rois', 'networks.nii')};

ANALYSIS_NAME = 'gPPI';

% 第一層只跑 DefaultMode 網路三個節點當 seed
ANALYSIS_SOURCES = { ...
    'networks.DefaultMode.MPFC', ...
    'networks.DefaultMode.PCC', ...
    'networks.DefaultMode.LP (L)' ...
};

BS_EFFECT_NAMES = {'AllSubjects', 'ExcludeOutlierSubjects'};
BS_CONTRAST     = [1 0];

CONTRAST.saveas                      = 'str_fst_r1_ran';
CONTRAST.between_conditions_names    = {'str_fst_r1', 'random'};
CONTRAST.between_conditions_contrast = [1 -1];

logFID = fopen(logFile, 'w');
write_log(logFID, '=== test_conn_save start ===\r\n\r\n');

%% 2. Setup + Denoising

batch = struct();
batch.filename     = batchFile;
batch.parallel.N   = 0;

batch.Setup.isnew                 = 1;
batch.Setup.nsubjects             = N_SUBJS;
batch.Setup.nsessions             = N_RUNS;
batch.Setup.RT                    = TR;
batch.Setup.conditions.names      = COND_NAMES;
batch.Setup.conditions.param      = zeros(1, N_CONDS);
batch.Setup.covariates.names      = COVAR_NAMES;
batch.Setup.subjects.effect_names = L2_COVARS;
batch.Setup.subjects.group_names  = {'AllSubjects'};
batch.Setup.subjects.groups       = ones(N_SUBJS, 1);
batch.Setup.subjects.add          = 1;

batch.Setup.rois.names = ROI_NAMES;
for iroi = 1:numel(ROI_NAMES)
    batch.Setup.rois.files{iroi} = ROI_FILES{iroi};
end
batch.Setup.rois.multiplelabels = 1;

meanMotionValues  = zeros(N_SUBJS, N_RUNS);
invalidScanCounts = zeros(N_SUBJS, N_RUNS);

for isub = 1:N_SUBJS
    sid  = SID_LIST(isub);
    subj = sprintf('sub-%03d', sid);

    batch.Setup.structurals{isub} = fullfile(bidsDir, subj, 'anat', ...
        sprintf('sub-%03d_space-%s_desc-preproc_T1w.nii.gz', sid, SPACE));
    batch.Setup.masks.Grey{isub}  = fullfile(bidsDir, subj, 'anat', ...
        sprintf('sub-%03d_space-%s_label-GM_probseg.nii.gz', sid, SPACE));
    batch.Setup.masks.White{isub} = fullfile(bidsDir, subj, 'anat', ...
        sprintf('sub-%03d_space-%s_label-WM_probseg.nii.gz', sid, SPACE));
    batch.Setup.masks.CSF{isub}   = fullfile(bidsDir, subj, 'anat', ...
        sprintf('sub-%03d_space-%s_label-CSF_probseg.nii.gz', sid, SPACE));

    for irun = 1:N_RUNS
        batch.Setup.functionals{isub}{irun} = fullfile(bidsDir, subj, 'func', ...
            sprintf('sub-%03d_task-srttprob_run-%02d_space-%s_desc-preproc_bold.nii.gz', sid, irun, SPACE));

        for icond = 1:N_CONDS
            condName = COND_NAMES{icond};
            [onsets, durations] = cond_setter(condName, irun, subj, bidsDir);
            batch.Setup.conditions.onsets{icond}{isub}{irun}    = onsets;
            batch.Setup.conditions.durations{icond}{isub}{irun} = durations;
        end

        for icov = 1:N_COVARS
            covarName = COVAR_NAMES{icov};
            covarFile = fullfile(bidsDir, subj, 'func', 'covariates', ...
                sprintf('%s_run-%02d.tsv', covarName, irun));
            covarData = readtable(covarFile, 'FileType', 'text', 'Delimiter', '\t');

            if ~isempty(covarData)
                batch.Setup.covariates.files{icov}{isub}{irun} = covarFile;
            end

            if matches(covarName, 'scrubbing')
                invalidScanCounts(isub, irun) = sum(covarData{:, 1});
            elseif matches(covarName, 'fd')
                meanMotionValues(isub, irun) = mean(covarData{:, 1});
            end
        end
    end
end

batch.Setup.subjects.effects{1} = mean(meanMotionValues, 2);
batch.Setup.subjects.effects{2} = sum(invalidScanCounts, 2);

batch.Setup.preprocessing.fwhm  = FWHM;
batch.Setup.preprocessing.steps = {'functional_smooth'};
batch.Setup.outputfiles = [1 1 0 0 0 0];
batch.Setup.done        = 1;
batch.Setup.overwrite   = 1;

batch.Denoising.confounds.names = CONFOUND_NAMES;
batch.Denoising.detrending      = POLY_ORD;
batch.Denoising.filter          = BP_HZ;

idx = find(strcmp(CONFOUND_NAMES, 'realignment'));
batch.Denoising.confounds.deriv{idx} = 1;
if MOT24, batch.Denoising.confounds.power{idx} = 2; end

idx = find(strcmp(CONFOUND_NAMES, 'aCompCor'));
batch.Denoising.confounds.dimensions{idx} = N_ACOMP;
idx = find(strcmp(CONFOUND_NAMES, 'aroma'));
batch.Denoising.confounds.dimensions{idx} = N_AROMA;

batch.Denoising.done      = 1;
batch.Denoising.overwrite = 1;

write_log(logFID, '[setup+denoise] start\r\n');
t0 = tic;
conn_batch(batch);
write_log(logFID, '[setup+denoise] done (%.1f sec)\r\n\r\n', toc(t0));

%% 3. QC scores —— 第一個重點測試

% 模擬修正後 run_conn_batch.m 的做法：用 conn('load',...) 取代 load + global
conn('load', batchFile);
global CONN_x;

write_log(logFID, '[QC] start\r\n');
t0 = tic;
s1 = conn_qascores('DataValidity',    [], []);
s2 = conn_qascores('DataQuality',     [], [], L2_COVARS, {});
s3 = conn_qascores('DataSensitivity', [], [], [], [], 'extreme');

% 關鍵補存：缺這行的話 ExcludeOutlierSubjects 只留在 memory
conn save;
write_log(logFID, '[QC] done (%.1f sec)  scores=[%.3f %.3f %.3f]\r\n\r\n', ...
    toc(t0), s1, s2, s3);

%% 4. First-level analysis —— 用 Analysis.sources 限定 3 顆 seed

write_log(logFID, '[1st-level] start  sources=%s\r\n', strjoin(ANALYSIS_SOURCES, ' | '));
t0 = tic;
conn_batch( ...
    'filename',             batchFile, ...
    'Analysis.name',        ANALYSIS_NAME, ...
    'Analysis.sources',     ANALYSIS_SOURCES, ...
    'Analysis.measure',     3, ...
    'Analysis.modulation',  1, ...
    'Analysis.conditions',  COND_NAMES, ...
    'Analysis.type',        3, ...
    'Analysis.done',        1, ...
    'Analysis.overwrite',   1 ...
);
write_log(logFID, '[1st-level] done (%.1f sec)\r\n\r\n', toc(t0));

%% 5. Second-level analysis —— 第二個重點測試

write_log(logFID, '[2nd-level] start: %s\r\n', CONTRAST.saveas);
t0 = tic;
conn_batch( ...
    'filename',                                 batchFile, ...
    'Results.saveas',                           CONTRAST.saveas, ...
    'Results.name',                             ANALYSIS_NAME, ...
    'Results.between_subjects.effect_names',    BS_EFFECT_NAMES, ...
    'Results.between_subjects.contrast',        BS_CONTRAST, ...
    'Results.between_conditions.effect_names',  CONTRAST.between_conditions_names, ...
    'Results.between_conditions.contrast',      CONTRAST.between_conditions_contrast, ...
    'Results.display',                          0, ...
    'Results.done',                             1 ...
);

% 關鍵補存：缺這行的話 Results.saved 只留在 memory
conn save;
write_log(logFID, '[2nd-level] done (%.1f sec)\r\n\r\n', toc(t0));

%% 6. 驗證：從 disk 重新載入，檢查兩個欄位是否真的存進去

clear CONN_x;
loaded = load(batchFile, 'CONN_x');
CONN_x = loaded.CONN_x;

l2names    = CONN_x.Setup.l2covariates.names;
savedNames = CONN_x.Results.saved.names;

hasExcl     = any(strcmp(l2names, 'ExcludeOutlierSubjects'));
hasContrast = any(strcmp(savedNames, CONTRAST.saveas));

write_log(logFID, '=== 驗證結果 ===\r\n');
write_log(logFID, 'Setup.l2covariates.names : %s\r\n', strjoin(l2names, ', '));
write_log(logFID, 'Results.saved.names      : %s\r\n', strjoin(savedNames, ', '));
write_log(logFID, '\r\n');
write_log(logFID, '[1] ExcludeOutlierSubjects 寫入 l2covariates : %s\r\n', tickmark(hasExcl));
write_log(logFID, '[2] %s 寫入 Results.saved              : %s\r\n', ...
    CONTRAST.saveas, tickmark(hasContrast));

if hasExcl && hasContrast
    write_log(logFID, '\r\nPASS\r\n');
else
    write_log(logFID, '\r\nFAIL\r\n');
end

fclose(logFID);

%% Helper functions

function [onsets, durations] = cond_setter(condName, irun, subj, bidsDir)
    onsets = []; durations = [];

    prefixes = {'str_', 'str_fst', 'str_snd', 'swi_'};
    if any(startsWith(condName, prefixes)) && ~contains(condName, num2str(irun))
        return;
    end

    if startsWith(condName, 'str_fst')
        fileTag = 'str_fst';
    elseif startsWith(condName, 'str_snd')
        fileTag = 'str_snd';
    elseif startsWith(condName, 'swi_')
        fileTag = 'switch';
    elseif startsWith(condName, 'str_')
        fileTag = 'structured';
    else
        fileTag = condName;
    end

    condFile = fullfile(bidsDir, subj, 'func', 'events', ...
        sprintf('%s_run-%02d.tsv', fileTag, irun));
    if isfile(condFile)
        condData  = readtable(condFile, 'FileType', 'text', 'Delimiter', '\t');
        onsets    = condData.onset;
        durations = condData.duration;
    end
end

function write_log(fid, fmt, varargin)
    msg = sprintf(fmt, varargin{:});
    fprintf('%s', msg);
    if fid > 0, fprintf(fid, '%s', msg); end
end

function s = tickmark(b)
    if b, s = 'OK'; else, s = 'MISSING'; end
end
