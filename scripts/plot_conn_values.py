#!/usr/bin/env python

import os
import argparse
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Config:
    def __init__(self, args):
        self.setup_paths(args)
        self.na_filler = 0
        self.seed = np.random.randint(0, 10000) if args.seed is None else args.seed
        self.cond_name_map = {
            "str_fst": "Sequence_1", 
            "str_snd": "Sequence_2", 
            "str"    : "Structured", 
            "swi"    : "Switch",
            "random" : "Random"
        }
        self.cond_colors = {
            "str_fst": "#006d77", 
            "str_snd": "#6a994e", 
            "str"    : "#006d77", 
            "swi"    : "#bc4749", 
            "random" : "#6c757d"
        }
        self.fig_dpi = 200
        self.fig_size = (10, 8)
        self.legend_loc = "lower center" # "upper center"
        self.legend_pos = (0.5, 1) # (0.5, -0.07)
        self.fig_w_ratios = [6, 1]
        self.title_fs = 18
        self.label_fs = 16
        self.ticks_fs = 14
        self.big_dot_size = 10
        self.small_dot_size = 16
        self.big_dot_alpha = 0.8
        self.small_dot_alpha = 0.5
        self.jitter_sd = 0.08
        self.errbar_capsize = 3
        self.errbar_lw = 1.5
        self.zero_lw = 0.6
        self.ctrl_lw = 1        
        
    def setup_paths(self, args):
        self.source_dir = os.path.dirname(os.path.abspath(__file__))
        self.xls_file = os.path.join(self.source_dir, "..", "conn_out", "values", args.project, f"ROI_{args.roi_prefix}.xlsx")
        self.fig_path_template_1 = os.path.join(self.source_dir, "..", "figures", "Connectivity", args.project, f"ROI_{args.roi_prefix}_{{ridx1}}_{{ridx2}}.png")
        self.fig_path_template_2 = os.path.join(self.source_dir, "..", "figures", "Connectivity", args.project, "mean_only", f"ROI_{args.roi_prefix}_{{ridx1}}_{{ridx2}}.png")
    

@dataclass
class CondData:
    name: str
    color: str
    runs: np.ndarray
    values: np.ndarray # shape: (n_subjects, n_runs)

    @property
    def mean(self) -> np.ndarray:
        return np.nanmean(self.values, axis=0)

    @property
    def se(self) -> np.ndarray:
        return np.nanstd(self.values, axis=0) / np.sqrt(np.sum(~np.isnan(self.values), axis=0))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot ROI-to-ROI connectivity strength across runs by condition.",
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        epilog="If -i (--roi_index_pairs) is not omitted, the script will print the mapping of ROI indices to ROI names and exit."
    )
    parser.add_argument("-p", "--project", type=str, required=True, 
                        help="Name of the CONN project.")
    parser.add_argument("-r", "--roi_prefix", type=str, default=["networks", "joanne"][1], 
                        help="Prefix of the ROI names used to generate the report.")
    parser.add_argument("-i", "--roi_index_pairs", type=str, nargs="+", default=None, 
                        help="One or more index pairs of the source and target ROIs. Format: 'source,target'. Example: -i 0,1 0,2 1,2")
    parser.add_argument("-s", "--seed", type=int, default=None, 
                        help="Random seed for plotting jitter.")
    parser.add_argument("-o", "--overwrite", action="store_true", 
                        help="Overwrite existing output files.")    
    return parser.parse_args()


def make_roi_dict(df: pd.DataFrame) -> dict[int, str]:
    '''
    Create a mapping of ROI indices to ROI names from the DataFrame.
    '''
    assert "roi2" in df.columns, "Expected 'roi2' column not found in the XLSX file."
    return {idx: name for idx, name in enumerate(df["roi2"].unique())}


def print_roi_dict(roi_dict: dict[int, str], roi_prefix: str):
    '''
    Print the mapping of ROI indices to ROI names.
    '''
    print(f"\nROI index for '{roi_prefix}':")
    print("\n".join([f"\t{idx:2d}:\t'{name}'" for idx, name in roi_dict.items()]))
    print("\nPlease specify the source and target ROI indices using the -i (--roi_index_pairs) arguments.")
    print()


def parse_ridx_pairs(raw_pairs: list[str]) -> list[tuple[int, int]]:
    '''
    Parse a list of strings representing ROI index pairs into a list of tuples.
    '''
    int_pairs = []
    for pair in raw_pairs:
        try:
            source, target = map(int, pair.split(","))
            int_pairs.append((source, target))
        except ValueError:
            raise ValueError(f"Invalid ROI index pair format: '{pair}'. Expected format: 'source,target' with integers.")
    return int_pairs


def parse_condition(df: pd.DataFrame, config: Config) -> pd.DataFrame:
    '''
    Parse the 'condition' column to extract the base condition name and run number, 
    and create new columns for them.
    '''
    assert "condition" in df.columns, "Expected 'condition' column not found in the XLSX file."
    extracted = df["condition"].str.extract(r"(.*)_r(\d+)$")
    df = df.assign(
        tmp_cond=extracted[0], 
        run=pd.to_numeric(extracted[1].fillna(config.na_filler), errors="coerce"),
        cond=lambda x: x["tmp_cond"].fillna(x["condition"])
    )
    return df.drop(columns=["tmp_cond", "condition"])


def make_cond_datas(df: pd.DataFrame, config: Config) -> dict[str, CondData]:
    '''
    Create a dictionary of CondData objects for each condition in the DataFrame.
    '''
    cond_datas = {}
    for cond_name, cond_df in df.groupby("cond"):
        wide_df = cond_df.pivot_table(index="subject", columns="run", values="value", aggfunc="mean")
        wide_df = wide_df.sort_index(axis=0).sort_index(axis=1)
        cond_datas[cond_name] = CondData(
            name=config.cond_name_map.get(cond_name, cond_name), 
            color=config.cond_colors.get(cond_name, "#000000"), 
            runs=wide_df.columns.to_numpy(), 
            values=wide_df.to_numpy()
        )
        
    return cond_datas


def plot_cond_data(ax: plt.Axes, cond_data: CondData, config: Config, plot_scatter: bool):
    '''
    Plot the connectivity values for a single condition on the given Axes.
    '''
    # Plot individual data points with jitter
    if plot_scatter:
        for run in cond_data.runs:
            y_arr = cond_data.values[:, cond_data.runs == run].flatten()
            rng = np.random.default_rng(seed=config.seed + run)
            x_jitter = rng.normal(0, config.jitter_sd, size=y_arr.size)
            ax.scatter(
                x=np.full_like(y_arr, fill_value=run) + x_jitter, 
                y=y_arr, 
                s=config.small_dot_size, 
                c=cond_data.color, 
                edgecolors="none", 
                alpha=config.small_dot_alpha, 
                zorder=1
            )
    
    # Plot mean and SE
    ax.errorbar(
        x=cond_data.runs, 
        y=cond_data.mean, 
        yerr=cond_data.se, 
        fmt="o-", 
        color=cond_data.color, 
        alpha=config.big_dot_alpha,
        lw=config.errbar_lw, 
        capsize=config.errbar_capsize, 
        markersize=config.big_dot_size,
        label=cond_data.name, 
        zorder=2
    )


def setup_axes(ax: plt.Axes, ax_id: int, run_arr: np.ndarray, ctrl_data: CondData, config: Config):
    '''
    Set up ticks, labels, spines, and reference lines for the given Axes based on its ID (left or right panel).
    '''
    if ax_id == 0: # for left panel
        ax.set_xlabel("Run", fontsize=config.label_fs)
        ax.set_ylabel("Connectivity Strength", fontsize=config.label_fs)
        ax.set_xticks(run_arr)
        ax.tick_params(axis="both", labelsize=config.ticks_fs)

    else: # for right panel
        ax.set_xlabel(f"\n{ctrl_data.name}", fontsize=config.label_fs)
        ax.set_xlim(config.na_filler - 0.5, config.na_filler + 0.5)
        ax.set_xticks([])
        ax.tick_params(axis="y", which="both", left=False, labelleft=False)
        ax.spines["left"].set_visible(False)

    # For both panels
    ax.axhline(0, color="black", ls="-", lw=config.zero_lw, zorder=3)
    if ctrl_data:
        ax.axhline(ctrl_data.mean, color=ctrl_data.color, ls="--", lw=config.ctrl_lw, zorder=3)


def plot_one_fig(cond_datas: dict[str, CondData], ctrl_cond: str, run_arr: np.ndarray, plot_scatter: bool, fig_title: str, fig_path: str, config: Config):
    '''
    Create and save a figure comparing connectivity values across conditions, 
    with an optional control condition in a separate panel.
    '''
    fig_cols = 2 if ctrl_cond else 1
    fig_kwargs = {"width_ratios": config.fig_w_ratios} if fig_cols == 2 else {}
    fig, axes = plt.subplots(1, fig_cols, sharey=True, **fig_kwargs)

    for cond_name, data in cond_datas.items():
        ax = axes if fig_cols == 1 else axes[1] if cond_name == ctrl_cond else axes[0]
        plot_cond_data(ax, data, config, plot_scatter=plot_scatter)
    
    # Add reference lines and finalize axes
    if ctrl_cond:
        for ax_id, ax in enumerate(axes):
            setup_axes(ax, ax_id, run_arr, cond_datas[ctrl_cond], config)
        ax = axes[0]
    else:
        setup_axes(axes, 0, run_arr, None, config)
        ax = axes
    
    ax.legend(
        loc=config.legend_loc, 
        bbox_to_anchor=config.legend_pos, # pin loc corner to the specific coordinate (x, y)
        ncol=len(cond_datas),
        fontsize=config.label_fs, 
        frameon=False
    )
    plt.suptitle(fig_title, fontsize=config.title_fs, fontweight="bold")
    plt.savefig(fig_path)
    plt.close(fig)
    print(f"Figure saved to: {fig_path}")


def run_one_iter(df: pd.DataFrame, roi_dict: dict, ridx_pair: tuple, config: Config, overwrite: bool):
    '''
    Create and save a figure for the specified pair of source and target ROIs.
    '''
    # Exit if the figure already exists and the overwrite flag is not enabled
    fig_path_1 = config.fig_path_template_1.format(ridx1=str(ridx_pair[0]), ridx2=str(ridx_pair[1]))
    fig_path_2 = config.fig_path_template_2.format(ridx1=str(ridx_pair[0]), ridx2=str(ridx_pair[1]))

    for fp in [fig_path_1, fig_path_2]:
        os.makedirs(os.path.dirname(fp), exist_ok=True)
        if os.path.exists(fp) and (not overwrite):
            print(f"\nFigure '{fp}' already exists and overwrite is not enabled. Skipping current iteration ...")
            return

    # Filter the DataFrame for the specified source and target ROIs
    roi1 = roi_dict[ridx_pair[0]]
    roi2 = roi_dict[ridx_pair[1]]
    targ_df = df.query(f"roi1 == '{roi1}' & roi2 == '{roi2}'")
    targ_df = targ_df.drop(columns=["roi1", "roi2"])
    fig_title = f"{roi1} x {roi2}"
    print(f"\nPlotting connectivity values for '{fig_title}' ...")
    
    # Create 'run' column from the 'condition' column
    targ_df = parse_condition(targ_df, config)
    run_arr = sorted(targ_df["run"].dropna().unique())
    if config.na_filler in run_arr:
        idx = run_arr.index(config.na_filler)
        run_arr = np.delete(run_arr, idx)

    # Create CondData objects for each condition
    cond_datas = make_cond_datas(targ_df, config)

    # Identify the control condition (if any)
    ctrl_cond = targ_df.query(f"run == {config.na_filler}")["cond"].unique()
    ctrl_cond = None if len(ctrl_cond) == 0 else ctrl_cond[0] if len(ctrl_cond) == 1 else ValueError(f"Multiple control conditions found: {ctrl_cond}")

    # Create and save the figures
    for fp, plot_scatter in zip([fig_path_1, fig_path_2], [True, False]):
        plot_one_fig(cond_datas, ctrl_cond, run_arr, plot_scatter, fig_title, fp, config)


def main():
    args = parse_args()
    config = Config(args)

    assert os.path.exists(config.xls_file), f"Excel file not found: {config.xls_file}"

    # If ROI index pairs are not specified, print the ROI index mapping and exit.
    if args.roi_index_pairs is None:
        mini_df = pd.read_excel(config.xls_file, usecols=["roi2"], nrows=10000)
        mini_df = mini_df.sort_values(by="roi2").reset_index(drop=True)
        roi_dict = make_roi_dict(mini_df)
        print_roi_dict(roi_dict, args.roi_prefix)
        return
    
    # Load the entire DataFrame and remove duplicate rows
    print(f"\nLoading data from: '{config.xls_file}' ...")
    df = pd.read_excel(config.xls_file)
    df = df.drop(columns=["contrast"])
    df = df.drop_duplicates().reset_index(drop=True)

    # Filter the DataFrame to include only rows corresponding to the specified source and target ROIs
    df = df.sort_values(by="roi2").reset_index(drop=True)
    roi_dict = make_roi_dict(df)

    # Set up plotting parameters at runtime
    plt.rcParams.update({
        "axes.spines.top": False, 
        "axes.spines.right": False, 
        "figure.constrained_layout.use": True, 
        "figure.figsize": config.fig_size, 
        "figure.dpi": config.fig_dpi
    })

    # Parse the specified ROI index pairs and create a figure for each pair
    ridx_pairs = parse_ridx_pairs(args.roi_index_pairs)
    print(f"\nParsed {len(ridx_pairs)} ROI index pair(s) to plot.")

    for ridx_pair in ridx_pairs:
        run_one_iter(df, roi_dict, ridx_pair, config, args.overwrite)

    print("\nAll done!\n")


if __name__ == '__main__':
    main()
