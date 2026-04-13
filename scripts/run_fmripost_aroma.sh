#!/bin/bash

PROJ_ROOT="/media/data3/Joanne_SRT_pw"
BIDS_DIR="/home/fmri2404/chiaoen/SrttNewV2/Nifti"
PREP_DIR="${PROJ_ROOT}/data/fmriprep"
POST_DIR="${PROJ_ROOT}/data/post_aroma"

WORK_DIR="${PROJ_ROOT}/work/post_aroma"
if [ ! -d $WORK_DIR ]; then mkdir -p $WORK_DIR; fi

for subj_dir in $BIDS_DIR/sub-*/; do
	subj=`basename $subj_dir`
	sid=`echo $subj | sed "s/sub-//g"`
	
	if [ ! -d $OUT_DIR/$subj ]; then
		# docker run -it --rm                                \
			# -v ${BIDS_DIR}:/data:ro                        \
			# -v /home/fmri2404/chiaoen/SrttNewV2/Nifti/derivatives/FreeSurfer:/fs_data:ro \
			# -v /usr/local/freesurfer/license.txt:/opt/freesurfer/license.txt \
			# -v ${PREP_DIR}:/out                            \
			# -v /home/fmri2404/chiaoen/SrttNewV2/work:/work \
			# nipreps/fmriprep:latest                        \
				# /data /out participant                     \
				# --participant-label $sid                   \
				# --skip_bids_validation                     \
				# --output-spaces T1w MNI152NLin2009cAsym MNI152NLin6Asym \
				# --force syn-sdc                            \
				# --fs-subjects-dir /fs_data                 \
				# --stop-on-first-crash                      \
				# --nthreads 10                              \
				# --omp-nthreads 10                          \
				# -w /work

		docker run -it --rm              \
			-v ${BIDS_DIR}:/data:ro      \
			-v ${PREP_DIR}:/prep         \
			-v ${POST_DIR}:/out          \
			-v ${WORK_DIR}:/work         \
			nipreps/fmripost-aroma:main  \
				/data /out participant   \
				--participant-label $sid \
				-t "srttprob"            \
				--skip_bids_validation   \
				-d fmriprep=/prep        \
				--nthreads 10            \
				--omp-nthreads 10        \
				-w /work
				
			# https://fmripost-aroma.readthedocs.io/latest/usage.html
	fi
done

