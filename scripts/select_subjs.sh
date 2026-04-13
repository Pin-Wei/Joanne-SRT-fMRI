#!/usr/bin/env bash

SRC_DIR="../data/fmriprep"
if [[ ! -d $SRC_DIR ]]; then
	mkdir -p $SRC_DIR
	cp -r "/home/fmri2404/chiaoen/SrttNewV2/Nifti/derivatives/fMRIPrep/sub-*" $SRC_DIR
fi

DST_DIR="../data/fmriprep_dropped"
if [[ ! -d $DST_DIR ]]; then mkdir -p $DST_DIR; fi

SUBJ_LIST_FILE="../data/meta/subj_list.txt"
SUBJ_LIST_PATH=$(dirname "$SUBJ_LIST_FILE")
if [[ ! -d $SUBJ_LIST_PATH ]]; then mkdir -p $SUBJ_LIST_PATH; fi
if [[ -f $SUBJ_LIST_FILE ]]; then rm $SUBJ_LIST_FILE; fi # overwrite if file exists

REMOVE_LIST=( 3 5 7 11 23 )
## 003, 007: Did not respond in more than 10% of the trials during the experiment, presumably because they fell asleep.
## 011: Dropped out of the experiment.
## 005, 023: Did not show learning of the first sequence in RTs.

for subj_dir in $SRC_DIR/sub-*/; do 
	subj=$(basename "${subj_dir}")
	sid=$(echo $subj | tr -dc '0-9' | awk 'sub("^0*", "")')
	add_to_list=true
	
	for sidx in ${REMOVE_LIST[@]}; do
		if [[ $sid -eq $sidx ]]; then
			mv $subj_dir $DST_DIR
			if [[ -f "$SRC_DIR/${subj}.html" ]]; then
				mv "$SRC_DIR/${subj}.html" $DST_DIR
			fi
			echo "Moved $subj to $DST_DIR"
			add_to_list=false
			break
		fi
	done
	
	if $add_to_list; then
		echo $subj >> $SUBJ_LIST_FILE
		echo "Added $subj to $SUBJ_LIST_FILE"
	fi
done