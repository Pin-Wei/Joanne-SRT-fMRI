#!/usr/bin/env bash
set -e

# This code is written for easy execution and logging
PYTHON="/home/aclexp/miniforge3/envs/pinwei/bin/python"

$PYTHON plot_conn_values.py \
	-p "no_PM_260524" -o -i 4,9 4,8 4,7 4,22 19,7
	
 # 0:	'L Amygdala'
 # 1:	'L Caudal Anterior Cingulate Cortex'
 # 2:	'L Caudate'
 # 3:	'L Cerebellum Cortex'
 # 4:	'L Granular Frontal Area (9)'
 # 5:	'L Inferior Parietal Lobule'
 # 6:	'L Insula'
 # 7:	'L Middle Frontal Area (46)'
 # 8:	'L Precentral Gyrus'
 # 9:	'L Primary Motor Cortex (4)'
# 10:	'L Putamen'
# 11:	'L Rostral Anterior Cingulate Cortex'
# 12:	'L Superior Parietal Lobule'
# 13:	'L Supplementary Motor Area (6)'
# 14:	'L Thalamus Proper'
# 15:	'R Amygdala'
# 16:	'R Caudal Anterior Cingulate Cortex'
# 17:	'R Caudate'
# 18:	'R Cerebellum Cortex'
# 19:	'R Granular Frontal Area (9)'
# 20:	'R Inferior Parietal Lobule'
# 21:	'R Insula'
# 22:	'R Middle Frontal Area (46)'
# 23:	'R Precentral Gyrus'
# 24:	'R Primary Motor Cortex (4)'
# 25:	'R Putamen'
# 26:	'R Rostral Anterior Cingulate Cortex'
# 27:	'R Superior Parietal Lobule'
# 28:	'R Supplementary Motor Area (6)'
# 29:	'R Thalamus Proper'

# 	-i 4,8 4,9 4,10 4,23 4,24 19,10 12,7 12,22 27,21 9,21 27,22 12,22 12,7

