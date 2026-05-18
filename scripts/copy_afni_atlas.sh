#!/usr/bin/env bash
set -e

# Copy AFNI-distributed atlases to local folder (../data/meta/)
# Usage: bash copy_afni_atlas.sh $IP
#
# 1. Freesurfer MNI2009c Desikan-Killiany (DK) parcellation
#	 - see also: /usr/local/freesurfer/FreeSurferColorLUT.txt
#
# 2. Brodmann-Pijnenburg atlas MNI 2009c 
#    - Pijnenburg, R., et al (2021). 
#      Myelo- and cytoarchitectonic microstructural and functional human cortical 
#      atlases reconstructed in common MRI space. 
#      NeuroImage, 239, 118274.
# 
# 3. JulichBrain 3.0 for MNI 2009c asymmetric space
#    - Amunts, K, Mohlberg, H, Bludau, S, & Zilles, K (2020).
#      Julich-Brain: A 3D probabilistic atlas of the human brain's cytoarchitecture.
#      Science, 369(6506), 988992.
#
# see all available atlases in file: AFNI_atlas_spaces.niml
	
IP=$1

ATLAS_FILES=( \
	"FS.afni.MNI2009c_asym.nii.gz" \
	"Brodmann_pijn_afni.nii.gz" \
	"Julich_MNI2009c.nii.gz" \
)
ATLAS_NAMES=( \
	"FS_atlas" \
	"Brodmann_atlas" \
	"Julich_atlas" \
)
TEMP_FILE="label_table.txt"

if [[ $# -eq 0 || $# -gt 2 ]]; then
    echo "Error: Should provide computer's IP number (23 or 37) to specify path to AFNI's binary files."
    echo "Usage: bash copy_fs_atlas.sh $IP"
	exit 1
else
	case "$1" in
		23)
			AFNI_DIR="/home/aclexp/abin"
			PYTHON="/home/aclexp/miniforge3/envs/pinwei/bin/python"
			;;
		37)
			AFNI_DIR="/usr/local/abin"
			PYTHON="/home/aclexp/mambaforge/envs/quanta/bin/python"
			;;
		*)
			echo "Error: path to AFNI's binary files for this IP number is not defined."
			exit 1
	        ;;
	esac
fi

DST_DIR="../data/meta"
if [ ! -d $DST_DIR ]; then mkdir -p $DST_DIR; fi

for i in ${!ATLAS_FILES[@]}; do
	orig=${ATLAS_FILES[$i]}
	nii=$DST_DIR/"${ATLAS_NAMES[$i]}.nii"
	nii_gz=$DST_DIR/"${ATLAS_NAMES[$i]}.nii.gz"
	txt=$DST_DIR/"${ATLAS_NAMES[$i]}.txt"
	
	if [ ! -f $nii ]; then
		if [ ! -f $nii_gz ]; then
			if [ ! -f $DST_DIR/$orig ]; then 
				echo "Coping ${AFNI_DIR}/${orig} ..."
				cp $AFNI_DIR/$orig $DST_DIR
			fi
			echo "Renaming ${DST_DIR}/${orig} to ${nii_gz} ..."
			mv $DST_DIR/$orig $nii_gz
		fi
		echo "Unzipping ${nii_gz} ..."
		gunzip $nii
	fi

	if [ ! -f $txt ]; then
		echo "Making label file ${txt} ..."
		3dinfo -labeltable $nii > $DST_DIR/$TEMP_FILE
		$PYTHON modify_label_table.py $DST_DIR/$TEMP_FILE $txt
		rm $DST_DIR/$TEMP_FILE
	fi
done

