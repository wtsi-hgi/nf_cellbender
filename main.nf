#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

VERSION = "0.0.2" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    cellbender__rb__get_input_cells;
    cellbender__remove_background;
    cellbender__remove_background__qc_plots;
    cellbender__remove_background__qc_plots_2;
    cellbender__gather_qc_input;
} from "./modules/core.nf"


// Set default parameters.
params.output_dir           = "nf-preprocessing"
params.help                 = false
// Default parameters for cellbender remove_background
params.cellbender_rb = [
    estimate_params_umis: [value: [
        expected_nemptydroplets_umi_cutoff: 0,
        method_estimate_ncells: 'dropletutils::barcoderanks::inflection',
        lower_bound_umis_estimate_ncells: 1000,
        method_estimate_nemptydroplets: 'dropletutils::barcoderanks::inflection',
        lower_bound_umis_estimate_nemptydroplets: 10,
        upper_bound_umis_estimate_nemptydroplets: 100,
        estimate_nemptydroplets_umi_subtract_factor: 25
    ]],
    epochs: [value: [200]],
    learning_rate: [value: [0.001, 0.0001]],
    fpr: [value: [0.01, 0.05]],
]
// Define the help messsage.
def help_message() {
    log.info """
    ============================================================================
     single cell preprocessing ~ v${VERSION}
    ============================================================================

    Runs basic single cell preprocessing

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Mandatory arguments:
        --file_paths_10x    Tab-delimited file containing experiment_id and
                            path_data_10xformat columns.

    Optional arguments:
        --output_dir        Directory name to save results to. (Defaults to
                            '${params.output_dir}')

        -params-file        YAML file containing analysis parameters. See
                            example in example_runtime_setup/params.yml.

    Profiles:
        local               local execution
        lsf                 lsf cluster execution
    """.stripIndent()
}


// Boot message - either help message or the parameters supplied.
if (params.help){
    help_message()
    exit 0
} else {
    log.info """
    ============================================================================
     single cell preprocessing ~ v${VERSION}
    ============================================================================
    file_paths_10x                : ${params.file_paths_10x}
    output_dir (output folder)    : ${params.output_dir}
    """.stripIndent()
    // A dictionary way to accomplish the text above.
    // def summary = [:]
    // summary['file_paths_10x'] = params.file_paths_10x
    // log.info summary.collect { k,v -> "${k.padRight(20)} : $v" }.join("\n")
}


// Initalize Channels.
// Channel: example init
// Channel
//     .fromPath( params.file_paths_10x )
//     .println()
// Channel: required files
// Channel
//     .fromPath(params.file_paths_10x)
//     .splitCsv(header: true, sep: "\t", by: 1)
//     .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format))}
//     .view()
channel__file_paths_10x = Channel
    .fromPath(params.file_paths_10x)
    .splitCsv(header: true, sep: "\t", by: 1)
    .map{row -> tuple(
        row.experiment_id,
        file("${row.data_path_10x_format}/barcodes.tsv.gz"),
        file("${row.data_path_10x_format}/features.tsv.gz"),
        file("${row.data_path_10x_format}/matrix.mtx.gz"),
        row.ncells_expected,
        row.ndroplets_include_cellbender
    )}
//n_experiments = file(params.file_paths_10x).countLines()


// Run the workflow.
workflow {
    main:

        // Prep the data for cellbender
        // This is silly, but to get the ouput directory structure that we
        // want, I need to pass all of the cellbender params here
        cellbender__rb__get_input_cells(
            params.output_dir,
            channel__file_paths_10x,
            params.cellbender_rb.estimate_params_umis.value
        )
        // Correct counts matrix to remove ambient RNA
        cellbender__remove_background(
            params.output_dir,
            cellbender__rb__get_input_cells.out.cb_input,
            params.cellbender_rb.epochs.value,
            params.cellbender_rb.learning_rate.value,
            params.cellbender_rb.zdim.value,
            params.cellbender_rb.zlayers.value,
            params.cellbender_rb.low_count_threshold.value,
            params.cellbender_rb.fpr.value
        )
        // Make some basic plots
        cellbender__remove_background__qc_plots(
            cellbender__remove_background.out.cb_plot_input
        )

    if (params.cellbender__remove_background__qc_plots_2.run_task) {
	log.info "params.cellbender__remove_background__qc_plots_2.run_task is set to true"
	
	if (params.cellbender__remove_background__qc_plots_2.file_paths_cellranger_filtered_10x) {
	    if (! file("${params.cellbender__remove_background__qc_plots_2.file_paths_cellranger_filtered_10x}").isEmpty()) {    
		log.info "params.cellbender__remove_background__qc_plots_2.file_paths_cellranger_filtered_10x is set to file path ${params.cellbender__remove_background__qc_plots_2.file_paths_cellranger_filtered_10x}"
		log.info "will run 'cellbender__remove_background__qc_plots_2' Nextflow task"
		
		// read in tab-delimited table
		//     params.file_paths_10x
		// columns 'experiment_id' and  'data_path_10x_format'
		// where 'data_path_10x_format' is raw cellranger dir .../raw_feature_bc_matrix
		ch_experimentid_paths10x_raw = Channel
		    .fromPath(params.file_paths_10x)
		    .splitCsv(header: true, sep: "\t", by: 1)
		    .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format))}

		// read in tab-delimited table
		//     params.cellbender__remove_background__qc_plots_2.file_paths_cellranger_filtered_10x
		// columns 'experiment_id' and  'data_path_10x_format'
		// where 'data_path_10x_format' is filtered cellranger dir .../filtered_feature_bc_matrix
		ch_experimentid_paths10x_filtered = Channel
		    .fromPath(params.cellbender__remove_background__qc_plots_2.file_paths_cellranger_filtered_10x)
		    .splitCsv(header: true, sep: "\t", by: 1)
		    .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format))}
		
		cellbender__remove_background.out.experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
		    .combine(ch_experimentid_paths10x_raw, by: 0)
		    .combine(ch_experimentid_paths10x_filtered, by: 0)
		    .combine(Channel.from("${params.cellbender_rb.fpr.value}"
					  .replaceFirst(/]$/,"")
					  .replaceFirst(/^\[/,"")
					  .split()))
		    .set{input_channel_qc_plots_2}
		
		cellbender__remove_background__qc_plots_2(input_channel_qc_plots_2)
		
	    } else {
		log.info "input parameter params.cellbender__remove_background__qc_plots_2.file_paths_cellranger_filtered_10x must be set to path to tsv file"
	    }
        }
    }
    
        // Bring all of the QC metrics together
        // Group by the second output which is the tuple
        // This code assumes tuple[0] = outdir
        //                   tuple[1] = experiment_id
        //                   tuple[2] = file to merge
        // gather_qc_input = cellbender__remove_background.out.filtered_10x.groupTuple(by: 0)
        //     .reduce([:]) { map, tuple ->
        //         def map_key = "epochs" + tuple[3][0]
        //         map_key = map_key + "__learnrt" + tuple[4][0].toString().replaceAll("-", "neg").replaceAll("\\.", "pt")
        //         map_key = map_key + "__zdim" + tuple[5][0]
        //         map_key = map_key + "__zlayer" + tuple[6][0]
        //         map_key = map_key + "__lowcount" + tuple[7][0]
        //         map_key = map_key + "__fpr" + tuple[8][0].toString().replaceAll("\\.", "pt").replaceAll(" ", "_")
        //         def key_list = map[map_key]
        //         if (!key_list) {
        //             key_list = [[tuple[1][0], tuple[2][0]]]
        //         } else {
        //             key_list.add([tuple[1][0], tuple[2][0]])
        //         }
        //         map[map_key] = key_list
        //         return(map)
        //     }
        //     .flatMap()
        //     .map {  entry ->
        //         combined_data = [entry.key, [], []]
        //         entry.value.each {
        //             combined_data[1].add(it[0]) // experiment_id
        //             combined_data[2].add(it[1]) //
        //         }
        //         return(combined_data)
        //     }
        cellbender__gather_qc_input(
            params.output_dir,
            cellbender__remove_background.out.results_list.collect()
        )
    // NOTE: One could do publishing in the workflow like so, however
    //       that will not allow one to build the directory structure
    //       depending on the input data call. Therefore, we use publishDir
    //       within a process.
    // publish:
    //     merge_samples.out.anndata to: "${params.output_dir}",
    //         mode: "copy",
    //         overwrite: "true"
}


workflow.onComplete {
    // executed after workflow finishes
    // ------------------------------------------------------------------------
    log.info """\n
    ----------------------------------------------------------------------------
     pipeline execution summary
    ----------------------------------------------------------------------------
    Completed         : ${workflow.complete}
    Duration          : ${workflow.duration}
    Success           : ${workflow.success}
    Work directory    : ${workflow.workDir}
    Exit status       : ${workflow.exitStatus}
    """.stripIndent()
}
