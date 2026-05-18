#!/usr/bin/env bash
set -e

FILE_PATH="/home/fmri2404/chiaoen/SrttNewV2/Nifti/code/SecondLevelAnalysis_ROI/ROI/ROI_masks"
ROI_FILES=$( ls $FILE_PATH | grep +tlrc )

DST_DIR="../data/roi_masks/by_joanne"
if [ ! -d $DST_DIR ]; then mkdir -p $DST_DIR; fi

for fn in ${ROI_FILES[@]}; do
	echo "Coping ${fn} ..."
	cp $FILE_PATH/$fn $DST_DIR
done
