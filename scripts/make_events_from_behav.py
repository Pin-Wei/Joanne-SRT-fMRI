#!/usr/bin/env python

import os
import glob
import shutil
import numpy as np
import pandas as pd

import warnings
warnings.simplefilter('ignore', np.RankWarning)

class Config:
    def __init__(self):
        self.setup_paths()
        self.sid_list = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25]
        self.run_list = range(1, 9)
        self.min_trials = 5
        self.seq_len = 8
        self.overwrite = False

    def setup_paths(self):
        self.source_dir = os.path.dirname(os.path.abspath(__file__))
        self.input_dir = os.path.join(self.source_dir, "..", "data", "behavioral")
        self.raw_data_path = os.path.join(self.input_dir, "raw_combined.csv")
        self.agg_data_path = os.path.join(self.input_dir, "rt_cleaned_medians.csv")
        self.out_data_path = os.path.join(self.input_dir, "rt_with_derivatives.csv")
        self.out_event_path_template = os.path.join(self.source_dir, "..", "data", "fmriprep", "{}", "func", "events", "{}_run-{}.tsv")

def compute_slopes(df, config):
    '''
    Compute slopes and intercepts for each structured block.
    Only correct trials are considered.
    If the number of correct trials is less than min_trials, the slope will be zero.
    If the slope is positive, the slope will be 0.
    If the intercept is negative, the intercept will be 0.
    For all random blocks, the slopes are 0, and the intercepts are the median RT.
    '''
    print(f"Participants and runs that have less than {config.min_trials} correct trials:")
    slope_dicts = []
    block = 0
    i = 0
    while i < df.shape[0]:
        sub_df = df.iloc[i:i+config.seq_len, :].reset_index(drop=True)
        assert len(set(sub_df["participant"])) == 1, "Should calc within the same participant"
        assert len(set(sub_df["session"])) == 1, "Should calc within the same run"
        assert len(set(sub_df["rule"])) == 1, "Should calc for the same rule"

        mask = sub_df["correct"].to_numpy()
        mask = np.where(mask == 0, False, True)
        x = np.arange(1, config.seq_len+1)[mask]
        y = sub_df["rt"].to_numpy()[mask]
        rule = sub_df["rule_group"][0]

        if mask.sum() < config.min_trials:
            print(f'sid: {sub_df["participant"][0]}, run: {sub_df["session"][0]}')
            slope, intercept = 0, 0
        elif rule == "Ran":
            slope, intercept = 0, np.median(y)
        else:
            slope, intercept = np.polyfit(x, y, 1)
            if slope > 0:
                slope = 0
                if intercept < 0:
                    intercept = 0
            
        slope_dicts.append({
            "participant": sub_df["participant"][0], 
            "session": sub_df["session"][0], 
            "block": block, 
            "rule_group": rule, 
            "slope": slope, 
            "intercept": intercept
        })
        i += config.seq_len
        block += 1
    
    return pd.DataFrame(slope_dicts)

def mark_str_switch(df, rule_col="rule_group", r1_name="Str1", r2_name="Str2"):
    '''
    Add columns to indicate:
    - "str_switch": whether the rule switched in the current trial
    - "str_switch_type": the type of the switch (1_to_2 or 2_to_1)
    '''
    df["prev_rule"] = df.groupby(["participant", "session"])[rule_col].shift(1)
    df["str_switch"] = (
        (df[rule_col] != df["prev_rule"]) & 
        df["prev_rule"].isin([r1_name, r2_name]) & 
        df[rule_col].isin([r1_name, r2_name])
    )
    df["str_switch_type"] = np.where(
        df["str_switch"] & (df["prev_rule"] == r1_name) & (df[rule_col] == r2_name), 
        "1_to_2", 
        np.where(
            df["str_switch"] & (df["prev_rule"] == r2_name) & (df[rule_col] == r1_name), 
            "2_to_1", 
            None
        )
    )
    return df

def compute_adapt_time(df):
    '''
    Compute the adapt time for the rule switched trials.
    The adapt time is the number of trials that the participant took to return to the baseline RT.
    If the participant took zero trials to return to the baseline RT, the adapt time will be the number of trials until the second time baseline RT is reached.
    '''
    def _compute_adapt_time_by_group(sub_df):
        sub_df = sub_df.reset_index()
        switch_idxs = sub_df.query("str_switch == True").index.tolist()
        for idx1, idx2 in zip(switch_idxs, switch_idxs[1:] + [-1]):
            rt_vals = sub_df.iloc[idx1:idx2]["rt"].values
            baseline = np.median(sub_df.iloc[idx1:idx2].query("correct == 1")["rt"].values)
            idxs = np.where(rt_vals <= baseline)[0]
            duration = idxs[0] if idxs[0] != 0 else idxs[1] 
            sub_df.loc[idx1, "adapt_time"] = duration
        return sub_df.set_index("index")

    df["adapt_time"] = np.nan
    df = (
        df.groupby(["participant", "session"])
        .apply(_compute_adapt_time_by_group, include_groups=False)
        .reset_index(level=[0, 1])
    )
    
    return df

def recreate_entire_folder(config):
    event_folder_template = os.path.dirname(config.out_event_path_template.format(f"sub-*", "*", "*"))
    event_folders = glob.glob(event_folder_template)

    if len(event_folders) > 0:
        print(f"\nFound {len(event_folders)} folders to remove.")
        event_folders = sorted(event_folders)
        
        # print(f"Do you want to recreate event folders? (y/n)\n{event_folders[0]}")
        # if input() not in ["y", "Y", "yes"]: 
        #     return        

    for folder in event_folders:
        print(f"Removing folder: {folder}")
        shutil.rmtree(folder, ignore_errors=True)
        os.makedirs(folder, exist_ok=True)

def make_event_files(df, config, default_dur=1, default_mod=1):
    for (sid, run), sub_df in df.groupby(["participant", "session"]):
        for k, trials in {
            "incorrect" : sub_df.query("correct == 0").reset_index(drop=True), 
            "random"    : sub_df.query("rule == 'Ran'").reset_index(drop=True), 
            "structured": sub_df.query("rule in ['StrA', 'StrB']").reset_index(drop=True), 
            "switch"    : sub_df.query("str_switch == 1").reset_index(drop=True), 
            # "adapt"     : sub_df.query("str_switch == 1").reset_index(drop=True), 
            "slope"     : sub_df, 
            "intercept" : sub_df
        }.items():
            fp = config.out_event_path_template.format(f"sub-{sid:03d}", k, f"{run:02d}")

            dur = {
                "switch": config.seq_len, 
                "adapt": trials["adapt_time"]
            }.get(k, default_dur)

            mod = {
                "slope": trials["slope"], 
                "intercept": trials["intercept"]
            }.get(k, default_mod)

            event_df = pd.DataFrame({
                "onset": trials["stim_onset"],
                "duration": dur,
                "modulation": mod
            }, index=trials.index)

            event_df.to_csv(fp, sep="\t", index=False, float_format="%.3f")
            print(f"File saved: {fp}")

# def copy_and_rename_files(config):
#     for event_name in ["structured", "switch"]:
#         for sid in config.sid_list:
#             for run in config.run_list:
#                 fp = config.out_event_path_template.format(f"sub-{sid:03d}", event_name, f"{run:02d}")
#                 run_group = {
#                     1: "r12", 2: "r12", 
#                     3: "r34", 4: "r34", 
#                     5: "r56", 6: "r56", 
#                     7: "r78", 8: "r78"
#                 }
#                 new_name = f"{event_name[:3]}-{run_group[run]}"
#                 fp2 = fp.replace(event_name, new_name)
#                 if not os.path.exists(fp2):
#                     print(f"Copying file to: {fp2}")
#                     shutil.copy(fp, fp2)

def main():
    config = Config()

    if config.overwrite:
        recreate_entire_folder(config)

    if not os.path.exists(config.out_data_path) or config.overwrite:
        data_raw = pd.read_csv(config.raw_data_path)
        data_raw = (
            data_raw.query(f"participant in {config.sid_list}")
            .reset_index(drop=True)
        )
        data_agg = pd.read_csv(config.agg_data_path)
        data = (
            data_raw.merge(
                data_agg, 
                on=["participant", "session", "cycle", "rule"],
                how="left"
            )
        )
        data = mark_str_switch(data)
        data = compute_adapt_time(data)
        data.insert(2, "block", np.repeat(np.arange(0, data.shape[0] // config.seq_len), config.seq_len))
        slope_df = compute_slopes(data, config)
        data = (
            data.merge(
                slope_df, 
                on=["participant", "session", "block", "rule_group"],
                how="left"
            )
        )
        data.to_csv(config.out_data_path, index=False)
    else:
        data = pd.read_csv(config.out_data_path)

    make_event_files(data, config)
    # copy_and_rename_files(config)

if __name__ == "__main__":
    main()

# ===================================================================================
# Participants and runs that have less than 5 correct trials:
# sid: 1, run: 5
# sid: 2, run: 6
# sid: 8, run: 1
# sid: 10, run: 4
# sid: 13, run: 1
# sid: 15, run: 3
# sid: 15, run: 7
# sid: 16, run: 6
# sid: 25, run: 1
# sid: 25, run: 1
# sid: 25, run: 1