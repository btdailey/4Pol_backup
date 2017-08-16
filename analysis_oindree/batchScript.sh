#!/bin/bash
cd /home/dailey.110/analysis_oindree

qsub -q long -l walltime=72:00:00 runvalgrind.sh
