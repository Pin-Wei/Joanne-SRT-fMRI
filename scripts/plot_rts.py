#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

class Config:
    def __init__(self):
        self.setup_paths()
        self.sid_list = [1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25]
        self.run_list = range(1, 9)

    def setup_paths(self):
        self.source_dir = os.path.dirname(os.path.abspath(__file__))
        self.input_dir = os.path.join(self.source_dir, "..", "data", "behavioral")
        self.aggdata_path = os.path.join(self.input_dir, "rt_cleaned_medians.csv")
        self.data_path = os.path.join(self.input_dir, "rt_with_derivatives.csv")
        self.out_fig_path_template = os.path.join(self.source_dir, "..", "figures", "{}", "sub-{}.png")

def plot_bars_with_color(
    ax, sub_data, x_col, y_col, 
    flip, minus_min, minus_mean, pred_y, 
    c_dict={0: "yellow", 1: "gray"}, c_col="correct", 
    bg_dict={"Ran": "white", "StrA": "red", "StrB": "green"}, bg_col="rule"
):
    '''
    Plot color-coded bars.
    The color of each bar is determined by the value of the 'c_col' column.
    (default: correct = gray, incorrect = yellow)
    The background color of each bar is determined by the value of the 'bg_col' column.
    (default: Ran = white, StrA = red, StrB = green)
    '''
    y = sub_data[y_col]
    if flip:
        y = 1 - y
    if minus_min:
        baseline = y.min()
    elif minus_mean:
        baseline = y.mean()
    else:
        baseline = 0
    y -= baseline
    x = np.arange(0, len(y)) if x_col == None else sub_data[x_col]

    c_list = []
    x0 = x[0]
    prev_label = sub_data[bg_col][0]
    
    for i, row in sub_data.iterrows():
        c = c_dict.get(row.get(c_col, 2), "tab:blue")
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
        
    return ax, baseline

def plot_agg_rt_bars(
    DF, config, folder_name, x_col=None, y_col="rt_medians", 
    flip=False, minus_min=False, minus_mean=False, 
    x_tick=10, y_tick=0.1, x_lim=(0, 256), y_lim=None, 
    label_fs=16, title_fs=18, fig_size=(20, 3), dpi=300, overwrite=False
):
    '''
    For each participant, plot the median RT for each run in one row.
    '''
    for sid in config.sid_list:
        fp = config.out_fig_path_template.format(folder_name, f"{sid:03d}")

        if os.path.exists(fp) and not overwrite:
            continue
        elif not os.path.exists(fp):
            os.makedirs(os.path.dirname(fp), exist_ok=True)

        plt.figure(figsize=fig_size, dpi=dpi)
        ax = plt.gca()
        sub_data = DF.query(f"participant == {sid}").reset_index(drop=True)
        ax, y_baseline = plot_bars_with_color(
            ax, sub_data, x_col, y_col, 
            flip, minus_min, minus_mean, pred_y=False
        )
        ax.xaxis.set_major_locator(ticker.MultipleLocator(x_tick))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(y_tick))
        ax.tick_params(axis="both", which="major", labelsize=label_fs)
        ax.set(xlim=x_lim, ylim=y_lim)
        ax.set_title(f"sub-{sid:03d}, baseline = {y_baseline:.3f}", fontsize=title_fs)
        plt.tight_layout()
        plt.savefig(fp, dpi=dpi)
        print(f"Figure saved: {fp}")
        plt.close()

def plot_rt_bars(
    DF, config, folder_name, x_col=None, y_col="rt", 
    flip=False, minus_min=False, minus_mean=False, pred_y=False, 
    x_tick=8, y_tick=0.1, x_lim=(0, 320), y_lim=None, 
    label_fs=16, title_fs=18, fig_size=(20, 20), dpi=300, overwrite=False
):
    '''
    For each participant, plot the RTs for each run in rows.
    Predicted RTs will be plotted as a line if 'pred_y' is True.
    '''
    for sid in config.sid_list:
        fp = config.out_fig_path_template.format(folder_name, f"{sid:03d}")

        if os.path.exists(fp) and not overwrite:
            continue
        elif not os.path.exists(fp):
            os.makedirs(os.path.dirname(fp), exist_ok=True)
        
        plt.figure(figsize=fig_size, dpi=dpi)

        for r, run in enumerate(config.run_list):
            ax = plt.subplot(len(config.run_list), 1, r+1)
            sub_data = DF.query(f"participant == {sid} & session == {run}").reset_index(drop=True)
            ax, y_baseline = plot_bars_with_color(
                ax, sub_data, x_col, y_col, 
                flip, minus_min, minus_mean, pred_y
            ) 
            ax.xaxis.set_major_locator(ticker.MultipleLocator(x_tick))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(y_tick))
            ax.tick_params(axis="both", which="major", labelsize=label_fs)
            ax.set(xlim=x_lim, ylim=y_lim)
            ax.set_title(f"run-{run:02d}, baseline = {y_baseline:.3f}", fontsize=title_fs)
            # ax.legend(bbox_to_anchor=(1.1, .5))
        
        plt.tight_layout()
        plt.savefig(fp, dpi=dpi)
        print(f"Figure saved: {fp}")
        plt.close()

def main():
    config = Config()

    agg_data_2 = pd.read_csv(config.aggdata_path)
    plot_agg_rt_bars(agg_data_2, config, folder_name="Median RTs (original)")
    plot_agg_rt_bars(agg_data_2, config, folder_name="Median RTs (flipped & minus min)", flip=True, minus_min=True)

    data = pd.read_csv(config.data_path)
    plot_rt_bars(data, config, folder_name="Raw RTs")
    plot_rt_bars(data, config, folder_name="Raw RTs (with slopes)", pred_y=True)

if __name__ == "__main__":
    main()