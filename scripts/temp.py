#!/usr/bin/env python

# This script copies specific covariate files (e.g., "fd" and "dvars") from the fmriprep output directory, renaming them with a "QC_" prefix. 
# It can also remove the copied files if needed.

import os
import glob
import shutil

def copy_specific_files(covar_name, new_name):
    files = glob.glob(os.path.join("..", "data", "fmriprep", "sub-*", "func", "covariates", f"{covar_name}_run-*.tsv"))

    for fp in files:
        fn = os.path.basename(fp)
        new_fn = fn.replace(covar_name, new_name)
        new_fp = os.path.join(os.path.dirname(fp), new_fn)
        print(f"Copying {fp} to {new_fp}")
        shutil.copy(fp, new_fp)

def remove_specific_files(covar_name):
    files = glob.glob(os.path.join("..", "data", "fmriprep", "sub-*", "func", "covariates", f"{covar_name}_run-*.tsv"))
    for fp in files:
        print(f"Removing file: {fp}")
        os.remove(fp)

if __name__ == "__main__":
    copy_specific_files("fd", "QC_fd")
    copy_specific_files("dvars", "QC_dvars")

    # remove_specific_files("QC_fd")
    # remove_specific_files("QC_dvars")