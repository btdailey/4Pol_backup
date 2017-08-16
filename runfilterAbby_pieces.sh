#!/bin/bash

#PBS -l walltime=34:00:00
#PBS -l nodes=1:fast:ppn=1
#PBS -j oe
#PBS -o /home/dailey.110/analysis_oindree/qsuboutputs/filter_study/

LOCAL_DIR=`pwd`
echo "Local dir is $LOCAL_DIR"

mkdir -p /tmp/pulser/$RUNNUMBER/$SIMFILE/$FILTER/$PHASE/$FILENAME

FILENAMER=/tmp/pulser/$RUNNUMBER/$SIMFILE/$FILTER/$PHASE/$FILENAME

cd /home/dailey.110/analysis_oindree
source /home/dailey.110/.bash_profile

./loopfulldata $RUNNUMBER $FILENAMER $FILTER $PHASE $start_number $end_number

if ((FILTER == 0)) 
then OUTPUT_DIR=no_fill
    elif((FILTER ==1))
then OUTPUT_DIR=Abby
    elif((FILTER == 2))
then OUTPUT_DIR=wiener
    elif((FILTER == 3))
then OUTPUT_DIR=geom_4pol_partial_0301
    elif((FILTER == 4)) 
then OUTPUT_DIR=no_filter
    elif((FILTER == 5))
then OUTPUT_DIR=subtraction
fi

if ((PHASE == 0))
then PHASE_DIR=full_run
    elif((PHASE == 1))
then PHASE_DIR=pulser_new
    elif((PHASE == 2))
then PHASE_DIR=pulser_response
    elif((PHASE == 3))
then PHASE_DIR=pulser
    elif((PHASE == 4))
then PHASE_DIR=thermal_shifted_mean_zero_simple_interp
fi
echo $OUTPUT_DIR
echo "exited and trying to move file!"

cp /tmp/pulser/$RUNNUMBER/$SIMFILE/$FILTER/$PHASE/$FILENAME/output* /data/anita/btdailey/final_filter/$PHASE_DIR/$OUTPUT_DIR
#mv $LOCAL_DIR/tmp/$CWDIR/$FILTER/$FILENAME/output* /data/anita/dailey.110/filter_study/saturated/noCW
rm -rf /tmp/pulser/$RUNNUMBER/$SIMFILE/$FILTER/$PHASE/


