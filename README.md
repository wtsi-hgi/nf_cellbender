
# Description

The methods used in this module are described in `docs/methods.pdf`. TODO: `docs/methods.pdf`

Below is the structure of the results directory. The values that will be listed in `description_of_params` within the directory structure correspond to the various parameters one can set. An example of a paramters file is found in `example_runtime_setup/params.yml`.
```bash
nf-qc_cluster
├── normalization_001::description_of_params
│   ├── [files: data]
│   ├── reduced_dims-pca::description_of_params
│   │   ├── [files: data]
│   │   ├── [plots: umap]
│   │   ├── cluster_001::description_of_params
│   │   │   ├── [files: data,clusters]
│   │   │   ├── [plots: umap]
│   │   │   ├── cluster_markers_001::description_of_params
│   │   │   │   ├── [files: cluster_marker_genes]
│   │   │   │   └── [plots: marker_genes,marker_genes_dotplot]
│   │   │   ├── cluster_markers_002::description_of_params
│   │   │   ... etc. ...
│   │   ├── cluster_002::description_of_params
│   │   ... etc. ...
│   ├── reduced_dims-harmony_001::description_of_params
│   ├── reduced_dims-harmony_002::description_of_params
│   ... etc. ...
├── normalization_002::description_of_norm_params
... etc. ...
└── adata.h5  # concatenated single cell data with no normalization
```


# TODO list

* Add `docs/methods.pdf` file.
* Add brief description of module.


# Enhancement list

* `scanpy_merge-dev.py`: If it were important to have a per sample filter, merge could be re-designed to accommodate this.
* `scanpy_cluster.py`: Currently for clustering, we can change method (leiden or louvain), resolution, and n_pcs. Are there other parameters that need to be scaled over?
* Check phenotypes against predicted sex from gene expression.
* Add basic QC plots - try to do this in R from anndata frame?
* Scrublet functionality + add to metadata + cluster distributions
* Gene scores + add to metadata
* Add marker gene AUC like here http://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters
* Add summary ARI and LISI metrics computed over a list of many different cluster annotations?
* Add tSNE plots - rapid plots with OpenTSNE?
* Calculate marker genes with diffxpy or logreg?


# Quickstart

Quickstart for deploying this pipeline locally and on a high performance compute cluster.


## 1. Set up the environment

Install the required packages via conda:
```bash
# The repo directory.
REPO_MODULE="${HOME}/repo/path/to/this/pipeline"

# Install environment using Conda.
conda env create --name sc_qc_cluster --file ${REPO_MODULE}/env/environment.yml

# Activate the new Conda environment.
source activate sc_qc_cluster

# To update environment file:
#conda env export --no-builds | grep -v prefix | grep -v name > environment.yml
```


## 2. Prepare the input files

Generate and/or edit input files for the pipeline.

The pipeline takes as input:
1. **--file_paths_10x**:  Tab-delimited file containing experiment_id and data_path_10x_format columns (i.e., list of input samples). Reqired.
2. **--file_metadata**:  Tab-delimited file containing sample metadata. This will automatically be subset down to the sample list from 1. Reqired.
3. **--file_sample_qc**:  YAML file containing sample qc and filtering parameters. Optional. NOTE: in the example config file, this is part of the YAML file for `-params-file`.
4. **--genes_exclude_hvg**:  Tab-delimited file with genes to exclude from
highly variable gene list. Must contain ensembl_gene_id column. Optional.
5. **--genes_score**:  Tab-delimited file with genes to use to score cells. Must contain ensembl_gene_id and score_idvcolumns. If one score_id == "cell_cycle", then requires a grouping_id column with "G2/M" and "S" (see example file in `example_runtime_setup`). Optional.
6. **-params-file**:  YAML file containing analysis parameters. Optional.
7. **--run_multiplet**:  Flag to run multiplet analysis. Optional.
8. **--file_cellmetadata**:  Tab-delimited file containing experiment_id and data_path_cellmetadata columns. For instance this file can be used to pass per cell doublet annotations. Optional.

Examples of all of these files can be found in `example_runtime_setup/`.


## 3. Set up and run Nextflow

Run Nexflow locally (NOTE: if running on a virtual machine you may need to set `export QT_QPA_PLATFORM="offscreen"` for scanpy as described [here](https://github.com/ipython/ipython/issues/10627)):
```bash
# Boot up tmux session.
tmux new -s nf

# Here we are not going to filter any variable genes, so don't pass a file.
# NOTE: All input file paths should be full paths.
nextflow run "${REPO_MODULE}/main.nf" \
    -profile "local" \
    --file_paths_10x "${REPO_MODULE}/example_runtime_setup/file_paths_10x.tsv" \
    --file_metadata "${REPO_MODULE}/example_runtime_setup/file_metadata.tsv" \
    --genes_score "${REPO_MODULE}/example_runtime_setup/genes_score_v001.tsv" \
    -params-file "${REPO_MODULE}/example_runtime_setup/params.yml"
```


Run Nextflow using LSF on a compute cluster. More on bgroups [here](https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsb.params.default_jobgroup.5.html).:
```bash
# Set the results directory.
RESULTS_DIR="/path/to/results/dir"
mkdir -p "${RESULTS_DIR}"

# Boot up tmux session.
tmux new -s nf

# Log into an interactive session.
# NOTE: Here we set the -G parameter due to our institute's LSF configuration.
bgadd "/${USER}/logins"
bsub -q normal -G team152 -g /${USER}/logins -Is -XF -M 8192 -R "select[mem>8192] rusage[mem=8192]" /bin/bash
# NOTE: If you are running over many cells, you may need to start an
# interactive session on a queue that allows long jobs
#bsub -q long -G team152 -g /${USER}/logins -Is -XF -M 18192 -R "select[mem>18192] rusage[mem=18192]" /bin/bash

# Activate the Conda environment (inherited by subsequent jobs).
conda activate sc_qc_cluster

# Set up a group to submit jobs to (export a default -g parameter).
bgadd -L 500 "/${USER}/nf"
export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
# Depending on LSF setup, you may want to export a default -G parameter.
export LSB_DEFAULTGROUP="team152"
# NOTE: By setting the above flags, all of the nextflow LSF jobs will have
# these flags set.

# Settings for scanpy (see note above).
export QT_QPA_PLATFORM="offscreen"

# Change to a temporary runtime dir on the node. In this demo, we will change
# to the same execution directory.
cd ${RESULTS_DIR}

# Remove old logs and nextflow output (if one previously ran nextflow in this
# dir).
rm -r *html;
rm .nextflow.log*;

# NOTE: If you want to resume a previous workflow, add -resume to the flag.
# NOTE: If you do not want to filter any variable genes, pass an empty file to
#       --genes_exclude_hvg. See previous local example.
# NOTE: --output_dir should be a full path - not relative.
nextflow run "${REPO_MODULE}/main.nf" \
    -profile "lsf" \
    --file_paths_10x "${REPO_MODULE}/example_runtime_setup/file_paths_10x.tsv" \
    --file_metadata "${REPO_MODULE}/example_runtime_setup/file_metadata.tsv" \
    --file_sample_qc "${REPO_MODULE}/example_runtime_setup/params.yml" \
    --genes_exclude_hvg "${REPO_MODULE}/example_runtime_setup/genes_remove_hvg_v001.tsv" \
    --genes_score "${REPO_MODULE}/example_runtime_setup/genes_score_v001.tsv" \
    --output_dir "${RESULTS_DIR}" \
    --run_multiplet \
    -params-file "${REPO_MODULE}/example_runtime_setup/params.yml" \
    -with-report "nf_report.html" \
    -resume

# NOTE: If you would like to see the ongoing processes, look at the log files.
cat .nextflow.log
```


Example of how one may sync results:
```bash
NF_OUT_DIR="/path/to/out_dir"
rsync -am --include="*.png" --include="*/" --exclude="*" my_cluster_ssh:${NF_OUT_DIR} .
rsync -am --include="*.png" --include="*/" --exclude="*" my_cluster_ssh:${NF_OUT_DIR} .
```

# Notes

* On 10 April 2020, we found nextflow was writing some output into the `${HOME}` directory and had used up the alotted ~15Gb on the Sanger farm. This resulted in a Java error as soon as a nextflow command was executed. Based on file sizes within `${HOME}`, it seemed like the ouput was being written within the conda environment (following `du -h | sort -V -k 1`). By deleting and re-installing the coda environment, the problem was solved. The below flags may help prevent this from the future. In addition, setting the flag `export JAVA_OPTIONS=-Djava.io.tmpdir=/path/with/enough/space/` may also help.

```bash
# To be run from the execution dir, before the above nextflow command
# If you are running this on a cluster, make sure you log into an interactive
# session with >25Gb of RAM.
export NXF_OPTS="-Xms25G -Xmx25G"
export NXF_HOME=$(pwd)
export NXF_WORK="${NXF_HOME}/.nexflow_work"
export NXF_TEMP="${NXF_HOME}/.nexflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/.nexflow_conda"
export NXF_SINGULARITY_CACHEDIR="${NXF_HOME}/.nexflow_singularity"
```
