#!/usr/bin/env python

import os
import glob
import shutil
import argparse
import numpy as np
import pandas as pd

from nilearn.glm.first_level import compute_regressor

class Config:
    def __init__(self):
        self.setup_paths()
        self.sid_list = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25]
        self.run_list = range(1, 9)
        self.tr = 2.0
        self.n_volumes = 200
        self.frame_times = np.arange(self.n_volumes) * self.tr

    def setup_paths(self):
        self.source_dir = os.path.dirname(os.path.abspath(__file__))
        self.prep_dir = os.path.join(self.source_dir, "..", "data", "fmriprep")
        self.event_path_template = os.path.join(self.prep_dir, "{}", "func", "events", "{}_run-{}.tsv")
        self.confound_path_template = os.path.join(self.prep_dir, "{}", "func", "{}_task-srttprob_run-{}_desc-confounds_timeseries.tsv")
        self.aroma_dir = os.path.join(self.source_dir, "..", "data", "post_aroma")
        self.aroma_path_template = os.path.join(self.aroma_dir, "{}", "func", "{}_task-srttprob_run-{}_desc-aroma_timeseries.tsv")
        self.out_path_template = os.path.join(self.prep_dir, "{}", "func", "covariates", "{}_run-{}.tsv")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--overwrite", action="store_true", 
                        help="Overwrite existing files.")
    parser.add_argument("-re", "--recreate", action="store_true", 
                        help="Recreate entire covariate folder for all subjects.")
    parser.add_argument("-rm", "--remove_files", type=str, nargs="*", 
                        help="Names of covariate files to remove.")
    return parser.parse_args()

def recreate_entire_folder(config):
    folder_template = os.path.dirname(config.out_path_template.format("sub-*", "*", "*"))
    folders = glob.glob(folder_template)

    if len(folders) > 0:
        print(f"\nFound {len(folders)} folders to remove.")
        folders = sorted(folders)
        
        # print(f"Do you want to recreate covariate folders? (y/n)\n{folders[0]}")
        # if input() not in ["y", "Y", "yes"]: 
        #     return        

    for folder in folders:
        print(f"Removing folder: {folder}")
        shutil.rmtree(folder, ignore_errors=True)
        os.makedirs(folder, exist_ok=True)

def remove_specific_files(covar_name, config):
    files = glob.glob(config.out_path_template.format("*", covar_name, "*"))
    if len(files) > 0:
        print(f"\nFound {len(files)} files to remove.")
        files = sorted(files)
        print(f"Are you sure you want to remove these files? (y/n)\n{files[0]}")
        if input() not in ["y", "Y", "yes"]: 
            return

    for fp in files:
        print(f"Removing file: {fp}")
        os.remove(fp)

def rename_specific_files(covar_name, new_name, config):
    files = glob.glob(config.out_path_template.format("*", covar_name, "*"))
    if len(files) > 0:
        print(f"\nFound {len(files)} files to rename.")
        files = sorted(files)
        print(f"Are you sure you want to rename these files? (y/n)\n{files[0]}")
        if input() not in ["y", "Y", "yes"]: 
            return

    for fp in files:
        fn = os.path.basename(fp)
        new_fn = fn.replace(covar_name, new_name)
        new_fp = os.path.join(os.path.dirname(fp), new_fn)
        print(f"Renaming {fp} to {new_fp}")
        os.rename(fp, new_fp)
        
def gen_from_event_files(sid, run, config, args, hrf_model="glover"):
    for k in ["slope", "intercept"]:
        out_fp = config.out_path_template.format(f"sub-{sid:03d}", k, f"{run:02d}")

        if not os.path.exists(out_fp) or args.overwrite:
            in_fp = config.event_path_template.format(f"sub-{sid:03d}", k, f"{run:02d}")
            df = pd.read_csv(in_fp, sep="\t")
            exp_condition = np.array([
                df["onset"].values, 
                df["duration"].values, 
                df["modulation"].values
            ])
            signal, _ = compute_regressor( # computes regressor sampled at frame_times
                exp_condition, 
                hrf_model, 
                config.frame_times
            )
            out_df = pd.DataFrame(signal)
            out_df.to_csv(out_fp, sep="\t", index=False, header=False)
            print(f"File saved: {out_fp}")

def extract_from_confound_timeseries(sid, run, config, args):
    '''
    Extract confound timeseries from fMRIPrep output and save them as separate files.
    - 'realignment': the 6 rigid-body head-motion parameters
    - 'Friston24': the 24 head-motion parameters proposed by Friston et al. (1996), including the 6 realignment parameters, their temporal derivatives, their squares, and the squares of their temporal derivatives
    - 'fd': framewise displacement (FD) 
    - 'dvars': standardized DVARS (D referring to temporal derivative of timecourses, VARS referring to RMS variance over voxels) 
    - 'scrubbing': motion outliers identified by fMRIPrep (binary flags)
    - 'aCompCor': the anatomical CompCor components (principal components of signals from white matter and CSF)
    - 'wm': the average signal within the cerebral white matter (WM) mask
    - 'csf': the average signal within the cerebro-spinal fluid (CSF) mask
    - 'global': the average signal within the whole-brain mask
    '''
    fp = config.confound_path_template.format(f"sub-{sid:03d}", f"sub-{sid:03d}", f"{run:02d}")
    confounds = pd.read_csv(fp, sep="\t")
    motion_cols = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]
    for k, cols in {
        "realignment": motion_cols, 
        "Friston24": motion_cols + [ f"{c}_derivative1" for c in motion_cols ] + [ f"{c}_power2" for c in motion_cols ] + [ f"{c}_derivative1_power2" for c in motion_cols ], 
        "fd": ["framewise_displacement"],
        "dvars": ["std_dvars"],
        "scrubbing": [ c for c in confounds.columns if c.startswith("motion_outlier")], 
        "aCompCor": [ c for c in confounds.columns if c.startswith("a_comp_cor") ],
        "wm": ["white_matter"], 
        "csf": ["csf"],
        "global": ["global_signal"] 
    }.items():
        out_fp = config.out_path_template.format(f"sub-{sid:03d}", k, f"{run:02d}")
        
        if not os.path.exists(out_fp) or args.overwrite:
            df = confounds.loc[:, cols]
            if k in ["fd", "dvars"]:
                df.iloc[0, :] = 0
            if k == "scrubbing":
                if len(cols) == 0:
                    df = pd.DataFrame(np.zeros((config.n_volumes, 1)))
                elif len(cols) > 1:
                    df = np.max(df, axis=1)
            # if df.shape[0] != config.n_volumes:
            #     raise ValueError(f"Number of volumes in {out_fp} should be {config.n_volumes}, but it is {df.shape[1]}")
            df.to_csv(out_fp, sep="\t", index=False, header=False)
            print(f"File saved: {out_fp}")

def extract_from_aroma_timeseries(sid, run, config, args):
    '''
    The confounds generated by fMRIPost-AROMA,
    which are ICA component time series classified as "rejected" by ICA-AROMA.
    Columns starting with aroma_orth_motion_ are the noise ICA component time series, after z-scoring and orthogonalization with respect to the signal ICA component time series.
    '''
    out_fp = config.out_path_template.format(f"sub-{sid:03d}", "aroma", f"{run:02d}")
    if not os.path.exists(out_fp) or args.overwrite:
        in_fp = config.aroma_path_template.format(f"sub-{sid:03d}", f"sub-{sid:03d}", f"{run:02d}")
        df = pd.read_csv(in_fp, sep="\t")
        df = df.loc[:, [ c for c in df.columns if c.startswith("aroma_orth_motion_") ]]
        df.to_csv(out_fp, sep="\t", index=False, header=False)
        print(f"File saved: {out_fp}")

def main():
    args = parse_args()
    config = Config()

    if args.recreate:
        recreate_entire_folder(config)

    if args.remove_files:
        for k in args.remove_files:
            remove_specific_files(k, config)

    for sid in config.sid_list:        
        for run in config.run_list:
            gen_from_event_files(sid, run, config, args)
            extract_from_confound_timeseries(sid, run, config, args)
            extract_from_aroma_timeseries(sid, run, config, args)

    print("\nDone.\n")

if __name__ == "__main__":
    main()