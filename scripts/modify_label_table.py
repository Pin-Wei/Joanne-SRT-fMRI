#!/usr/bin/env python

import sys
import pandas as pd

input_path = sys.argv[1]
output_path = sys.argv[2]

df = pd.read_csv(input_path, sep="\\s+", header=None, index_col=0, skiprows=3)
df = df[:-1]
df.index = df.index.astype(int)
df = df.sort_index()
df.to_csv(output_path, sep="\t", header=False)

## Explore the image file ===========================================================

# import nibabel as nib
# from nilearn.image import math_img
# from nilearn.plotting import plot_roi, show

# img_path = "/usr/local/abin/FS.afni.MNI2009c_asym.nii.gz"
# img = nib.load(img_path)
# img_data = img.get_fdata()

# print("Labels: ", np.unique(img_data))

# roi_index = 6
# roi_name = df.loc[roi_index, 1]
# img_masked = math_img(f"img == {roi_index}", img=img)
# plot_roi(img_masked, title=f"ROI {roi_index}: {roi_name}")
# show()