#!/usr/bin/python

import os
import sys
import numpy as np
import pandas as pd
# import nilearn as nl
from nilearn.glm.first_level import compute_regressor

class Config:
    def __init__(self):
        self.sid_list = [1, 2, 4]
        self.run_list = range(1, 9)
        self.tr = 2.0
        self.n_volumes = 200
        self.source_dir = os.path.dirname(os.path.abspath(__file__))
        self.event_path_template = os.path.join(self.source_dir, "..", "data", "fmriprep_mini", "{}", "func", "events", "{}_run-{}.tsv")
        self.confound_path_template = os.path.join(self.source_dir, "..", "data", "fmriprep_mini", "{}", "func", "{}_task-srttprob_run-{}_desc-confounds_timeseries.tsv")
        self.out_path_template = os.path.join(self.source_dir, "..", "data", "fmriprep_mini", "{}", "func", "covariates", "{}_run-{}.tsv")
        self.overwrite = True # False

def extract_from_confound_timeseries(sid, run, config):
    fp = config.confound_path_template.format(f"sub-{sid:03d}", f"sub-{sid:03d}", f"{run:02d}")
    confounds = pd.read_csv(fp, sep="\t")
    for k, cols in {
        "motion": ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"], 
        "aCompCor": [ c for c in confounds.columns if c.startswith("a_comp_cor") ],
        "global": ["csf", "white_matter"] # "global_signal"
    }.items():
        out_fp = config.out_path_template.format(f"sub-{sid:03d}", k, f"{run:02d}")
        if not os.path.exists(out_fp) or config.overwrite:
            os.makedirs(os.path.dirname(out_fp), exist_ok=True)
            df = confounds.loc[:, cols]
            df.to_csv(out_fp, sep="\t", index=False)
            print(f"File saved: {out_fp}")

def gen_from_event_files(sid, run, config, hrf_model="glover"):
    frame_times = np.arange(config.n_volumes) * config.tr
    out_fp = config.out_path_template.format(f"sub-{sid:03d}", "behav", f"{run:02d}")
    if not os.path.exists(out_fp) or config.overwrite:
        os.makedirs(os.path.dirname(out_fp), exist_ok=True)
        out_df = pd.DataFrame()
        for k in ["slope", "intercept"]:
            fp = config.event_path_template.format(f"sub-{sid:03d}", k, f"{run:02d}")
            df = pd.read_csv(fp, sep="\t")
            exp_condition = np.array([
                df["onset"].values, 
                df["duration"].values, 
                df["modulation"].values
            ])
            signal, _ = compute_regressor( # computed regressor sampled at frame_times
                exp_condition, 
                hrf_model, 
                frame_times
            )
            out_df[k] = signal.squeeze()
        out_df.to_csv(out_fp, sep="\t", index=False)
        print(f"File saved: {out_fp}")

def main():
    config = Config()
    for sid in config.sid_list:
        for run in config.run_list:
            extract_from_confound_timeseries(sid, run, config)
            gen_from_event_files(sid, run, config)

if __name__ == "__main__":
    main()