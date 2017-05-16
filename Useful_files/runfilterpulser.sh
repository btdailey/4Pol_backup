#!/bin/bash

#PBS -l walltime=10:00:00
#PBS -l nodes=1:fast:ppn=2

#PBS -m abe
#PBS -M dailey.110@osu.edu
#PBS -j oe
#PBS -o /home/dailey.110/analysis_oindree/qsuboutputs/filter_study/
LOCAL_DIR=`pwd`
echo "Local dir is $LOCAL_DIR"

mkdir -p /tmp/$FILTER/$FILENAME

FILENAMER=/tmp/$FILTER/$FILENAME

cd /home/dailey.110/analysis_oindree
source /home/dailey.110/.bash_profile

./looppulserdata $RUNNUMBER $FILENAMER $FILTER $PHASE

if ((FILTER == 0)) 
then OUTPUT_DIR=no_fill
elif((FILTER ==1))
then OUTPUT_DIR=rayleigh
elif((FILTER == 2))
then OUTPUT_DIR=wiener
elif((FILTER == 3))
then OUTPUT_DIR=geom_4pol
elif((FILTER == 4)) 
then OUTPUT_DIR=no_filter
elif((FILTER == 5))
then OUTPUT_DIR=subtraction
fi

if ((PHASE == 0))
then PHASE_DIR=no_phase_change
    elif((PHASE == 1))
then PHASE_DIR=random_phase
    elif((PHASE == 2))
then PHASE_DIR=interpolated_phase
    elif((PHASE == 3))
then PHASE_DIR=geometric_method
    elif((PHASE == 4))
then PHASE_DIR=shifted_method ##NOT IMPLEMENTED IN CURRENT CODE
fi

echo $OUTPUT_DIR
echo "exited and trying to move file!"

cp /tmp/$FILTER/$FILENAME/output* /data/anita/btdailey/final_filter/$PHASE_DIR/$OUTPUT_DIR
#mv $LOCAL_DIR/tmp/$CWDIR/$FILTER/$FILENAME/output* /data/anita/btdailey/filter_study/saturated/noCW
rm -rf /tmp/$FILTER/$FILENAME


