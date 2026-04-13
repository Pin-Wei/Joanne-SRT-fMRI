#!/usr/bin/env bash

DST_DIR="../data/meta"
NEW_NAME="FS.afni_atlas"
TEMP_FILE="label_table.txt"

AFNI_DIR="/usr/local/abin"
ATLAS_FILE="FS.afni.MNI2009c_asym.nii.gz"

# AFNI-distributed Freesurfer DK parcellation atlas for MNI 2009c asym template
# see also: /usr/local/freesurfer/"FreeSurferColorLUT.txt

if [ ! -f $DST_DIR/"${NEW_NAME}.nii.gz" ]; then
	echo "Coping atlas image file ..."
	if [ ! -f $DST_DIR/$ATLAS_FILE ]; then cp $AFNI_DIR/$ATLAS_FILE $DST_DIR; fi
	mv $DST_DIR/$ATLAS_FILE $DST_DIR/"${NEW_NAME}.nii.gz"
fi

if [ ! -f $DST_DIR/"${NEW_NAME}.txt" ]; then
	echo "Making atlas label file ..."
	3dinfo -labeltable $DST_DIR/"${NEW_NAME}.nii.gz" > $DST_DIR/$TEMP_FILE
	python modify_label_table.py $DST_DIR/$TEMP_FILE $DST_DIR/"${NEW_NAME}.txt"
	rm $DST_DIR/$TEMP_FILE
fi

