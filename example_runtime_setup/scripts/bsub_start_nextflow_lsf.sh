#!/usr/bin/env bash

# clean up previous run files
rm -f *.log
rm -f bsub.o
rm -f bsub.e
rm -f bjob.id

# start Nextflow via bsub:
bsub -G hgi \
     -R'select[mem>8000] rusage[mem=8000] span[hosts=1]' \
     -M 8000 \
     -n 2 \
     -o bsub.o -e bsub.e \
     -q basement \
     bash scripts/start_nextflow_lsf.sh > bjob.id

# get process PID 
echo "Nextflow Bjob ID saved in file bjob.id" 
echo kill with \"bkill ID_number\" command
echo "check logs files bsub.o, bsub.e and .nextflow.log"
