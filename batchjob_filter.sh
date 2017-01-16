#!/bin/bash
cd ~/analysis_oindree/

#qsub -q long -l walltime=72:00:00 batchfilterpulser.sh
#qsub -q long -l walltime=72:00:00 batchfilterthermal.sh
#qsub -q long -l walltime=72:00:00 batch10sample.sh
#qsub -q long -l walltime=72:00:00 batchsimulationdata.sh


qsub -q long -l walltime=72:00:00 batchsimulationdata_pieces.sh
#qsub -q long -l walltime=72:00:00 batchfilterpulser_pieces.sh
#qsub -q long -l walltime=72:00:00 batch10sample_pieces.sh
#qsub -q long -l walltime=72:00:00 batchfilterthermal_pieces.sh
#qsub -q long -l walltime=72:00:00 batch90sample_pieces.sh