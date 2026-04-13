#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from nilearn.glm.first_level import compute_regressor

class Config:
    def __init__(self):
        self.setup_paths()
        self.sid_list = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25]
        self.run_list = range(1, 9)
        self.tr = 2
        self.n_volumes = 200
        self.frame_times = np.arange(self.n_volumes) * self.tr
        self.overwrite = False

    def setup_paths(self):
        self.source_dir = os.path.dirname(os.path.abspath(__file__))
        self.input_path_template = os.path.join(self.source_dir, "..", "data", "fmriprep", "sub-{}", "func", "{}", "{}_run-{}.tsv")
        self.out_fig_path_template = os.path.join(self.source_dir, "..", "figures", "{}", "sub-{}.png")

def event_to_signls(
    event_df, config, hrf_model="glover"
):
    exp_condition = np.array([
        event_df["onset"].values, 
        event_df["duration"].values, 
        event_df["modulation"].values
    ])
    signal, _ = compute_regressor(
        exp_condition, 
        hrf_model, 
        config.frame_times
    )

    return signal

def make_event_timeseries(event_df, config):
    '''
    Create a time series data at units of TRs
    based on the onsets, durations, and modulations of the event.
    '''
    stim_times = np.zeros(config.tr * config.n_volumes) # unit in seconds
    for _, row in event_df.iterrows():
        t0 = row["onset"] - 1 # 0-indexed
        t1 = t0 + row.get("duration", 1)
        stim_times[t0:t1+1] = row.get("modulation", 1)
    stim_times = stim_times[::config.tr] # unit in TRs

    return stim_times

def plot_onsets(
    folder_name, data_names, data_types, config, 
    x_tick=20, label_fs=16, title_fs=16,
    fig_size=(15, 15), dpi=300, overwrite=False
):
    for sid in config.sid_list:
        out_fp = config.out_fig_path_template.format(folder_name, f"{sid:03d}")

        if os.path.exists(out_fp) and not overwrite:
            continue
        elif not os.path.exists(out_fp):
            os.makedirs(os.path.dirname(out_fp), exist_ok=True)

        plt.figure(figsize=fig_size)

        for r, run in enumerate(config.run_list):
            ax = plt.subplot(len(config.run_list), 1, r+1)
            x_pos_list, y_pos_list = [], []
            y_name_list = []
            y_pos = 0

            for data_name, data_type in zip(data_names, data_types):
                in_fp = config.input_path_template.format(f"{sid:03d}", data_type, data_name, f"{run:02d}")
                if data_type == "events":
                    event_df = pd.read_csv(in_fp, sep="\t")
                    timepoints = event_df["onset"].values
                elif data_type == "covariates":
                    stim_times = pd.read_csv(in_fp, sep="\t", header=None)[0].values
                    timepoints = np.where(stim_times != 0)[0]
                    timepoints = np.multiply(timepoints, config.tr)
                else:
                    raise ValueError(f"Unknown data type: {data_type}")
                
                count = len(timepoints)
                x_pos_list.append(timepoints)
                y_pos_list.append(y_pos)
                y_name_list.append(f"{data_name} ({count})")
                y_pos += 1

            ax.eventplot(x_pos_list, lineoffsets=y_pos_list)
            ax.set_yticks(y_pos_list, y_name_list)
            ax.set_xlim(0, config.tr * config.n_volumes)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(x_tick))
            ax.tick_params(axis="both", which="major", labelsize=label_fs)
            ax.set_title(f"run-{run:02d}", fontsize=title_fs)
        
        plt.tight_layout()
        plt.savefig(out_fp, dpi=dpi)
        print(f"Figure saved: {out_fp}")
        plt.close()

def plot_signals(
    folder_name, data_name, data_type, config, 
    label_fs=16, title_fs=16,
    fig_size=(15, 15), dpi=300, overwrite=False
):
    for sid in config.sid_list:
        out_fp = config.out_fig_path_template.format(folder_name, f"{sid:03d}")

        if os.path.exists(out_fp) and not overwrite:
            continue
        elif not os.path.exists(out_fp):
            os.makedirs(os.path.dirname(out_fp), exist_ok=True)

        plt.figure(figsize=fig_size)

        for r, run in enumerate(config.run_list):
            ax = plt.subplot(len(config.run_list), 1, r+1)

            in_fp = config.input_path_template.format(f"{sid:03d}", data_type, data_name, f"{run:02d}")
            if data_type == "events":
                event_df = pd.read_csv(in_fp, sep="\t")
                stim_times = make_event_timeseries(event_df, config)
                signal = event_to_signls(event_df, config)
            elif data_type == "covariates":
                signal = pd.read_csv(in_fp, sep="\t", header=None)[0].values
            else:
                raise ValueError(f"Unknown data type: {data_type}")
            
            ax.axhline(0, color="red", linestyle="--")
            if data_type == "events":
                ax.plot(config.frame_times, stim_times)
            ax.plot(config.frame_times, signal, color="black")

            ax.set_xlim(0, config.tr * config.n_volumes)
            ax.tick_params(axis="both", which="major", labelsize=label_fs)
            ax.set_title(f"run-{run:02d}", fontsize=title_fs)

        plt.tight_layout()
        plt.savefig(out_fp, dpi=dpi)
        print(f"Figure saved: {out_fp}")
        plt.close()

def main():
    config = Config()
    plot_onsets(
        folder_name="Bad timepoints", 
        data_names=["incorrect", "scrubbing"], 
        data_types=["events", "covariates"], 
        config=config
    )
    plot_signals(
        folder_name="Switching signals", 
        data_name="switch", 
        data_type="events", 
        config=config
    )
    plot_signals(
        folder_name="Structured block signals", 
        data_name="structured", 
        data_type="events", 
        config=config
    )
    plot_signals(
        folder_name="Random block signals", 
        data_name="random", 
        data_type="events", 
        config=config
    )
    print("\nDone!\n")

if __name__ == "__main__":
    main()