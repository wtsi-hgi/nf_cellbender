#!/usr/bin/env bash

# activate Nextflow conda env
conda init bash
eval "$(conda shell.bash hook)"
conda activate nextflow

# run nextflow main.nf with inputs and lsf config:
export NXF_OPTS="-Xms5G -Xmx5G"
nextflow run ./nf_cellbender/main.nf \
      --file_sample_qc params.yml -params-file params.yml \
      -c ./nf_cellbender/nextflow.config -c inputs.nf -profile lsf_hgi \
      --nf_ci_loc $PWD -resume
