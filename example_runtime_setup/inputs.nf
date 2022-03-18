params {
    
    output_dir = "${projectDir}/../outputs"
    
    // experiment_id, data_path_10x_format, data_path_raw_h5
    file_paths_10x = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/wbc_mult_donor/results/raw.Submission_Data_Pilot_UKB.file_paths_10x.tsv"

    cellbender__remove_background__qc_plots_2 {
	run_task = true // whether to run 'cellbender__remove_background__qc_plots_2 Nextflow' task.
	// provide Cellranger filtered output h5 (to be compared with cellbender output)
	// tab-delimited file, with 2 columns required in header (first row): experiment_id, data_path_filt_h5
	file_paths_cellranger_filtered_10x = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/wbc_mult_donor/results/Submission_Data_Pilot_UKB.file_paths_10x.tsv"
    }
}

    //file_paths_10x = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/fetch_gsheet_sbw_test/results/raw.Submission_Data_Pilot_UKB.file_paths_10x.tsv"
    //file_paths_10x = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/fetch_gsheet_sbw_test/results/cellbender_input.tsv"
    //file_paths_10x = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/fetch_gsheet_sbw_test/results/Submission_Data_Pilot_UKB.file_paths_10x.tsv"
//    output_dir           = "nf-qc_cluster"
