#!/usr/bin/env bash
set -e

ATLAS_LABELS=( DK BA )
ATLAS_DIR="../data/meta"
ROIMASK_DIR="../data/roi_masks"

for atlas_label in ${ATLAS_LABELS[@]}; do
	case $atlas_label in
		DK)
			atlas_nii="${ATLAS_DIR}/FS_atlas.nii"
			roi_indices=( \
				  5   6   7   8   9  15  48  53  69  71  74  80 \
				 25  26  27  28  29  32  83  88 104 106 109 115 \
			)
			roi_names=( \
				DK.L-Cerebellum-WM     \
				DK.L-Cerebellum-Cortex \
				DK.L-Thalamus          \
				DK.L-Caudate           \
				DK.L-Putamen           \
				DK.L-Amygdala          \
				DK.L-Caudal-ACC        \
				DK.L-Inferior-Parietal \
				DK.L-Precentral        \
				DK.L-Rostral-ACC       \
				DK.L-Superior-Parietal \
				DK.L-Insula            \
				DK.R-Cerebellum-WM     \
				DK.R-Cerebellum-Cortex \
				DK.R-Thalamus          \
				DK.R-Caudate           \
				DK.R-Putamen           \
				DK.R-Amygdala          \
				DK.R-Caudal-ACC        \
				DK.R-Inferior-Parietal \
				DK.R-Precentral        \
				DK.R-Rostral-ACC       \
				DK.R-Superior-Parietal \
				DK.R-Insula            \
			)
			;;
		BA)
			atlas_nii="${ATLAS_DIR}/Brodmann_atlas.nii"
			roi_indices=( \
				  3   5   8  38 \
				103 105 108 138 \
			)
			roi_names=( \
				BA.L-4  \
				BA.L-6  \
				BA.L-9  \
				BA.L-46 \
				BA.R-4  \
				BA.R-6  \
				BA.R-9  \
				BA.R-46 \
			)
			;;
		*)
			echo "Error: '${atlas}' not defined."
			exit 1
	        ;;
	esac
	
	for i in ${!roi_indices[@]}; do
		roimask_nii=$ROIMASK_DIR/${roi_names[$i]}.nii
		
		if [ ! -f $roimask_nii ]; then
			echo "==================================="
			echo "Making $(basename $roimask_nii) ..."
	
			mri_binarize \
				--i $atlas_nii --match ${roi_indices[$i]} \
				--o $roimask_nii
		fi
	done
done