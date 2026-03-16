#!/usr/bin/python

import os
import numpy as np
import pandas as pd
# import shutil
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import warnings
warnings.simplefilter('ignore', np.RankWarning)

class Config:
    def __init__(self):
        self.source_dir = os.path.dirname(os.path.abspath(__file__))
        self.input_dir = os.path.join(self.source_dir, "..", "data", "behavioral")
        self.raw_data_path = os.path.join(self.input_dir, "raw_combined.csv")
        self.agg_data_path = os.path.join(self.input_dir, "rt_cleaned_medians.csv")
        self.out_event_path_template = os.path.join(self.source_dir, "..", "data", "fmriprep_mini", "{}", "func", "events", "{}_run-{}.tsv")
        self.out_fig_path_template = os.path.join(self.source_dir, "..", "figures", "sub-{}.png")
        self.sid_list = [1, 2, 4]
        # self.sid_list = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25]

def compute_slopes(df, min_trials=5):
    '''
    Compute slopes and intercepts for each structured block.
    The length of each block is 8.
    Only correct trials are considered.
    If the number of correct trials is less than min_trials, the slope will be zero.
    If the slope is positive, the slope will be 0.
    If the intercept is negative, the intercept will be 0.
    For all random blocks, the slopes are 0, and the intercepts are the median RT.
    '''
    print(f"Participants and runs that have less than {min_trials} correct trials:")
    slope_dicts = []
    block = 0
    i = 0
    while i < df.shape[0]:
        sub_df = df.iloc[i:i+8, :].reset_index(drop=True)
        assert len(set(sub_df["participant"])) == 1, "Should calc within the same participant"
        assert len(set(sub_df["session"])) == 1, "Should calc within the same run"
        assert len(set(sub_df["rule"])) == 1, "Should calc for the same rule"

        mask = sub_df["correct"].to_numpy()
        mask = np.where(mask == 0, False, True)
        x = np.arange(1, 9)[mask]
        y = sub_df["rt"].to_numpy()[mask]
        rule = sub_df["rule_group"][0]

        if mask.sum() < min_trials:
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
        i += 8
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

def make_event_files(df, out_path):
    for (sid, run), sub_df in df.groupby(["participant", "session"]):
        # folder = os.path.dirname(out_path.format(f"sub-{sid:03d}", "*", "*"))
        # shutil.rmtree(folder, ignore_errors=True)

        for k, trials in {
            "incorrect" : sub_df.query("correct == 0").reset_index(drop=True), 
            "random"    : sub_df.query("rule == 'Ran'").reset_index(drop=True), 
            "structured": sub_df.query("rule in ['StrA', 'StrB']").reset_index(drop=True), 
            "switch"    : sub_df.query("str_switch == 1").reset_index(drop=True), 
            "adapt"     : sub_df.query("str_switch == 1").reset_index(drop=True), 
            "slope"     : sub_df, 
            "intercept" : sub_df
        }.items():
            fp = out_path.format(f"sub-{sid:03d}", k, f"{run:02d}")
            os.makedirs(os.path.dirname(fp), exist_ok=True)
            
            dur = {
                "adaptation": trials["adapt_time"]
            }.get(k, 1)

            mod = {
                "slope": trials["slope"], 
                "intercept": trials["intercept"]
            }.get(k, 1)

            event_df = pd.DataFrame({
                "onset": trials["stim_onset"],
                "duration": dur,
                "modulation": mod
            }, index=trials.index)

            event_df.to_csv(fp, sep="\t", index=False, float_format="%.3f")
            print(f"File saved: {fp}")

def plot_bars(
    DF, sid_list, run_list, y_col="rt", x_col=None, pred_y=False, 
    c_dict={0: "yellow", 1: "gray"}, c_col="correct", 
    bg_dict={"Ran": "white", "StrA": "red", "StrB": "green"}, bg_col="rule", 
    x_tick=8, x_lim=(0, 320), y_lim=None, 
    fig_size=(15, 15), dpi=300, out_path=None
):
    '''
    For each participant, plot the RTs for each run.
    The color of each bar is determined by the value of the c_col column.
    The background color of each bar is determined by the value of the bg_col column.
    Predicted RTs will be plotted as a line if pred_y is True.
    '''
    for sid in sid_list:
        plt.figure(figsize=fig_size, dpi=dpi)

        for r, run in enumerate(run_list):
            ax = plt.subplot(len(run_list), 1, r+1)

            sub_data = DF.query(f"participant == {sid} & session == {run}").reset_index(drop=True)
            y = sub_data[y_col]
            x = np.arange(0, len(y)) if x_col == None else sub_data[x_col]

            c_list = []
            x0 = x[0]
            prev_label = sub_data[bg_col][0]
            
            for i, row in sub_data.iterrows():
                c = c_dict.get(row[c_col], "tab:blue")
                c_list.append(c)
                
                curr_label = row[bg_col]
                if curr_label != prev_label:
                    x1 = x[i] -.5
                    color = bg_dict[prev_label]
                    ax.axvspan(x0, x1, facecolor=color, alpha=.5, label=prev_label)
                    x0 = x1
                prev_label = curr_label

            ax.bar(x, y, width=.7, color=c_list)
            
            if pred_y:
                sub_data["seq_order"] = np.tile(np.arange(1, 9), sub_data.shape[0] // 8)
                sub_data["rt_pred"] = sub_data["seq_order"] * sub_data["slope"] + sub_data["intercept"]
                
                for b in np.arange(x[0], x[-1], 8):
                    x_ = np.arange(b, b+8)
                    y_ = sub_data.iloc[b:b+8]["rt_pred"]
                    ax.plot(x_, y_, color="k")
                
            ax.xaxis.set_major_locator(ticker.MultipleLocator(x_tick))
            ax.set(xlim=x_lim, ylim=y_lim, title=f"run-{run:02d}")
            # ax.legend(bbox_to_anchor=(1.1, .5))
        
        plt.tight_layout()
        fp = out_path.format(f"{sid:03d}")
        os.makedirs(os.path.dirname(fp), exist_ok=True)
        plt.savefig(fp, dpi=dpi)
        print(f"Figure saved: {fp}")
        plt.close()

def main():
    config = Config()
    data_raw = pd.read_csv(config.raw_data_path)
    data_raw = data_raw.query(f"participant in {config.sid_list}").reset_index(drop=True)

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
    data.insert(2, "block", np.repeat(np.arange(0, data.shape[0] // 8), 8))
    slope_df = compute_slopes(data)
    data = (
        data.merge(
            slope_df, 
            on=["participant", "session", "block", "rule_group"],
            how="left"
        )
    )
    make_event_files(data, config.out_event_path_template)
    # plot_bars(
    #     data, config.sid_list, range(1, 9), pred_y=True, 
    #     out_path=config.out_fig_path_template
    # )

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