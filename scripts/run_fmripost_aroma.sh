#!/bin/bash

PROJ_ROOT="/media/data3/Joanne_SRT_pw"
BIDS_DIR="/home/fmri2404/chiaoen/SrttNewV2/Nifti"
PREP_DIR="${PROJ_ROOT}/data/fmriprep"
POST_DIR="${PROJ_ROOT}/data/post_aroma"

WORK_DIR="${PROJ_ROOT}/work/post_aroma"
if [ ! -d $WORK_DIR ]; then mkdir -p $WORK_DIR; fi

## argument count check --------------------------------------------------------------------

if [ $# -eq 0 ]; then
    echo "Error: Should provide at least one argument to specify which program to run."
    echo "(0: fmriprep, 1: fmripost-aroma, 2: both)."
    exit 1

elif [ $# -gt 2 ]; then
    echo "Error: Too many arguments provided."
    echo "argv[1]: 0 (fmriprep), 1 (fmripost-aroma), or 2 (both)."
    echo "argv[2]: int (subject ID)."
    exit 1

elif [ $# -eq 2 ]; then
    sid=$(printf "%03d" "$2")
    SUBJ_DIRS=( $BIDS_DIR/sub-${sid} )

else
    SUBJ_DIRS=( $BIDS_DIR/sub-* )
fi

## specify which program(s) to execute -----------------------------------------------------

case "$1" in
    0)
        RUN_PREP=true
        RUN_POST=false
        ;;
    1)
        RUN_PREP=false
        RUN_POST=true
        ;;
    2)
        RUN_PREP=true
        RUN_POST=true
        ;;
    *)
        echo "Error: The argument should be 0 (fmriprep), 1 (fmripost-aroma), or 2 (both)."
        exit 1
        ;;
esac

## main ------------------------------------------------------------------------------------

for subj_dir in ${SUBJ_DIRS[*]}; do
	subj=`basename $subj_dir`
	sid=`echo $subj | sed "s/sub-//g"`
	
	if [[ $RUN_PREP && ! -d $PREP_DIR/$subj ]]; then
		docker run -it --rm                                \
			-v ${BIDS_DIR}:/data:ro                        \
			-v /home/fmri2404/chiaoen/SrttNewV2/Nifti/derivatives/FreeSurfer:/fs_data:ro \
			-v /usr/local/freesurfer/license.txt:/opt/freesurfer/license.txt \
			-v ${PREP_DIR}:/out                            \
			-v /home/fmri2404/chiaoen/SrttNewV2/work:/work \
			nipreps/fmriprep:latest                        \
				/data /out participant                     \
				--participant-label $sid                   \
				--skip_bids_validation                     \
				--output-spaces T1w MNI152NLin2009cAsym MNI152NLin6Asym \
				--force syn-sdc                            \
				--fs-subjects-dir /fs_data                 \
				--stop-on-first-crash                      \
				--nthreads 10                              \
				--omp-nthreads 10                          \
				-w /work
	fi
	
	if [[ $RUN_POST && ! -d $OUT_DIR/$subj ]]; then
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
	fi
done

