#!/usr/bin/env nextflow

// NOTE: label 'big_mem' may be useful at some point


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process cellbender__rb__get_input_cells {
    // Calculates thresholds for input cells of cellbender__remove_background
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.equalsIgnoreCase("barcodes.tsv.gz")) {
                        null
                    } else if(filename.equalsIgnoreCase("features.tsv.gz")) {
                        null
                    } else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            val(ncells_expected),
            val(ndroplets_include_cellbender)
        )
        val(estimate_params_umis)

    output:
        val(outdir, emit: outdir)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            path("${outfile}-expected_cells.txt"),
            path("${outfile}-total_droplets_included.txt"),
            emit: cb_input
        )
        path(
            "${outfile}-expected_cells.txt",
            emit: expected_cells
        )
        path(
            "${outfile}-total_droplets_included.txt",
            emit: total_droplets_include
        )
        path("${outfile}-cell_estimate_cutoff.tsv.gz")
        path("${outfile}-total_droplets_cutoff.tsv.gz")
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/${experiment_id}"
        outdir = "${outdir}/cellbender-estimate_ncells_nemptydroplets"
        outfile = "umi_count_estimates"
        cmd__expected_ncells = ""
        if ("${ncells_expected}" != "NA") {
            cmd__expected_ncells = "--expected_ncells ${ncells_expected}"
            outdir = "${outdir}-ncells_expected__${ncells_expected}"
        }
        cmd__droplets_include = ""
        if ("${ndroplets_include_cellbender}" != "NA") {
            cmd__droplets_include = "--expected_nemptydroplets ${ndroplets_include_cellbender}"
            outdir = "${outdir}-nemptydrops_expected__${ndroplets_include_cellbender}"
        }
        if ("${outdir}" == "${outdir_prev}/${experiment_id}/cellbender-params_setbyuser") {
            outdir = "${outdir_prev}/${experiment_id}/cellbender-estimated_ncells_nemptydroplets"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellbender__remove_background__get_input_cells: ${process_info}"
        echo "outdir: ${outdir}"
        rm -fr plots
        mkdir txd_input
        ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
        ln --physical ${file_10x_features} txd_input/features.tsv.gz
        ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz

        python ${projectDir}/bin/015-get_estimates_from_umi_counts.py \\
            --tenxdata_path txd_input \\
            --output_file ${outfile} \\
            --expected_nemptydroplets_umi_cutoff \\
                ${estimate_params_umis.expected_nemptydroplets_umi_cutoff} \\
            --method_estimate_ncells \\
                ${estimate_params_umis.method_estimate_ncells} \\
            --lower_bound_umis_estimate_ncells \\
                ${estimate_params_umis.lower_bound_umis_estimate_ncells} \\
            --method_estimate_nemptydroplets \\
                ${estimate_params_umis.method_estimate_nemptydroplets} \\
            --lower_bound_umis_estimate_nemptydroplets \\
                ${estimate_params_umis.lower_bound_umis_estimate_nemptydroplets} \\
            --upper_bound_umis_estimate_nemptydroplets \\
                ${estimate_params_umis.upper_bound_umis_estimate_nemptydroplets} \\
            --estimate_nemptydroplets_add_umifactor \\
                ${estimate_params_umis.estimate_nemptydroplets_umi_add_factor} \\
            --estimate_nemptydroplets_subtract_dropletfactor \\
                ${estimate_params_umis.estimate_nemptydroplets_subtract_cell_factor} \\
            --estimate_nemptydroplets_min_nemptydroplets \\
                ${estimate_params_umis.estimate_nemptydroplets_min_drop} \\
                ${cmd__expected_ncells} ${cmd__droplets_include}

        # Clean up plots
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process cellbender__remove_background {
    // Remove ambient RNA
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    //maxForks 2         // hard to control memory usage. limit to 3 concurrent
    label 'gpu'          // use GPU
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename == "barcodes.tsv.gz") {
                        null
                    } else if(filename.equalsIgnoreCase("features.tsv.gz")) {
                        null
                    } else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            path(expected_cells),
            path(total_droplets_include)
        )
        each epochs
        each learning_rate
        each zdims
        each zlayers
        each low_count_threshold
        each fpr

    output:
        val(outdir, emit: outdir)
        tuple(
            val(outdir),
            val(experiment_id),
            path("${outfile}-filtered_10x_mtx-file_list.tsv"),
            val(epochs),
            val(learning_rate),
            val(zdims),
            val(zlayers),
            val(low_count_threshold),
            val(fpr),
            emit: cb_results_details
        )
        tuple(
            val(outdir),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            path("*_filtered.h5"),
            emit: cb_plot_input
        )
        path(
            "${runid}-${outfile}-filtered_10x_mtx-file_list.tsv",
            emit: results_list
        )
        path("${outfile}.log", emit: log)
        path("${outfile}_cell_barcodes.csv", emit: barcodes)
        path("*_unfiltered.h5", emit: unfiltered_h5s)
        path("*_filtered.h5", emit: filtered_h5s)
        path("*filtered_10x_mtx/barcodes.tsv.gz", emit: tenx_barcodes)
        path("*filtered_10x_mtx/features.tsv.gz", emit: tenx_features)
        path("*filtered_10x_mtx/matrix.mtx.gz", emit: tenx_matrix)
        path(file_10x_barcodes, emit: raw_tenx_barcodes)
        path(file_10x_features, emit: raw_tenx_features)
        path(file_10x_matrix, emit: raw_tenx_matrix)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true
        tuple(
            val(experiment_id),
            val(outdir),
            path("*_unfiltered.h5"),
            path(expected_cells),
            path(total_droplets_include),
            emit: experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
        )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/${experiment_id}"
        lr_string = "${learning_rate}".replaceAll("\\.", "pt")
        lr_string = "${lr_string}".replaceAll("-", "neg")
        fpr_string = "${fpr}".replaceAll("\\.", "pt").replaceAll(" ", "_")
        cb_params = "cellbender_params"
        cb_params = "${cb_params}-epochs_${epochs}"
        cb_params = "${cb_params}__learnrt_${lr_string}"
        cb_params = "${cb_params}__zdim_${zdims}"
        cb_params = "${cb_params}__zlayer_${zlayers}"
        cb_params = "${cb_params}__lowcount_${low_count_threshold}"
        outdir = "${outdir}/${cb_params}".replaceAll(
            "cellbender_params",
            "cellbender"
        )
        outfile = "cellbender"
        // outdir = "${outdir}-fpr_${fpr_string}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellbender__remove_background: ${process_info}"
        echo "outdir: ${outdir}"
        rm -fr plots
        mkdir txd_input
        ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
        ln --physical ${file_10x_features} txd_input/features.tsv.gz
        ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz
        # Run cellbender
        cellbender remove-background \\
            --input txd_input \\
            --output ${outfile} \\
            --cuda \\
            --expected-cells \$(cat ${expected_cells}) \\
            --total-droplets-included \$(cat ${total_droplets_include}) \\
            --model full \\
            --z-dim ${zdims} \\
            --z-layers ${zlayers} \\
            --low-count-threshold ${low_count_threshold} \\
            --epochs ${epochs} \\
            --learning-rate ${learning_rate} \\
            --fpr ${fpr}
        # If outfile does not have h5 appended to it, move it.
        [ -f ${outfile} ] && mv ${outfile} ${outfile}.h5
        # Clean up file names (e.g., FPR_0.1 becomes FPR_0pt1)
        # *.h5 will grep unfiltered and filtered
        # NOTE: We assume the only h5 files in the dir are from cellbender
        for i in \$(ls *.h5); do
            echo \$i
            out_file=\$(echo \$i | sed s/"\\."/"pt"/g | sed s/"pth5"/"\\.h5"/)
            # If outfile does not exist, move i to out_file
            [ ! -f \$out_file ] && mv \$i \$out_file
            echo \$out_file
        done
        # Make a table of cellbender output matricies in reference to their
        # final outdir... this will be used to make a group of files
        # for input
        032-clean_cellbender_results.py \\
            --nf_outdir_tag ${outdir} \\
            --cb_outfile_tag ${outfile} \\
            --experiment_id ${experiment_id} \\
            --fpr '${fpr}' \\
            --cb_params ${cb_params}
        # Move the list of output matricies to one with a unique job id so
        # that we can collect them all later in the pipeline across all job
        # iterations and not overwrite anything.
        cp ${outfile}-filtered_10x_mtx-file_list.tsv \\
            ${runid}-${outfile}-filtered_10x_mtx-file_list.tsv
        # Move plots to a plot dir
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process cellbender__remove_background__qc_plots {
    // QC plots from cellbdender
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    //maxForks 2         // hard to control memory usage. limit to 3 concurrent
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.equalsIgnoreCase("barcodes.tsv.gz")) {
                        null
                    } else if(filename.equalsIgnoreCase("features.tsv.gz")) {
                        null
                    } else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        tuple(
            val(outdir_prev),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            file(h5_filtered_cellbender)
        )

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}" // /${experiment_id}"
        h5_filtered_cellbender = h5_filtered_cellbender.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellbender__remove_background__qc_plots: ${process_info}"
        echo "outdir: ${outdir}"
        mkdir -p txd_input
        ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
        ln --physical ${file_10x_features} txd_input/features.tsv.gz
        ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz
        # Make a file with list of our files
        echo "${h5_filtered_cellbender}" | sed s/","/"\\n"/g > files.txt
        for i in \$(cat files.txt); do
            echo \$i
            out_file=\$(echo \$i | sed s/".h5"//)
            035-analyse_cellbender_results.py \\
                --tenxdata_path txd_input \\
                --h5_cellbender \$i \\
                --output_file cellbender_results-\$out_file \\
                --number_cpu ${task.cpus}
        done
        rm files.txt
        mkdir -p plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process cellbender__remove_background__qc_plots_2 {
    // Second set of QC plots from cellbdender.
    // This task compare the cellbender output with both the cellranger filtered and cellragner raw outputs
    // ------------------------------------------------------------------------
    tag { "$experiment_id" }
    publishDir "${params.output_dir}/$experiment_id/compare_cellranger/", pattern: "fpr_${fpr}/${experiment_id}/*.png", 
        saveAs: {filename ->
        filename.replaceAll("fpr_${fpr}/${experiment_id}/", "fpr_${fpr}/")
    },
        mode: "${task.publish_mode}",
        overwrite: "true"
    
    //cache false        // cache results from run
    //maxForks 2         // hard to control memory usage. limit to 3 concurrent
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script
    input:
        tuple val(experiment_id), val(outdir), path(cellbender_unfiltered_h5s), path(expectedcells), path(totaldropletsinclude), path(raw_cellranger_mtx), path(filtered_cellranger_mtx), val(fpr)
    output:
        val(outdir, emit: outdir)
        path("fpr_${fpr}/${experiment_id}/*.png"), emit: plots_png 

    script:
"""
fprid=\$(echo $fpr | sed s'/\\./pt/'g)
cellbender_unfiltered_h5=cellbender_FPR_\${fprid}_unfiltered.h5

echo fpr $fpr
echo fprid \$fprid
echo cellbender_unfiltered_h5 \$cellbender_unfiltered_h5

n_expected_cells=\$(cat $expectedcells)
echo n_expected_cells \${n_expected_cells}
n_total_droplets_included=\$(cat $totaldropletsinclude)
echo n_total_droplets_included \${n_total_droplets_included}

python ${projectDir}/bin/037-plot_cellranger_vs_cellbender.py \\
    --samplename \"${experiment_id}\" \\
    --raw_cellranger_mtx \"${raw_cellranger_mtx}\" \\
    --filtered_cellranger_mtx \"${filtered_cellranger_mtx}\" \\
    --cellbender_unfiltered_h5 \"\$cellbender_unfiltered_h5\" \\
    --fpr \"${fpr}\" \\
    --n_expected_cells \"\${n_expected_cells}\" \\
    --n_total_droplets_included \"\${n_total_droplets_included}\" \\
    --out_dir \$PWD
"""
}


process cellbender__gather_qc_input {
    // Prepare cell bender output for qc_cluster pipeline. For each epoch and
    // learning rate, and fpr, gather 10x matrix files into format for
    // qc_cluster.
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        file(cb_results_tsvs)

    output:
        val(outdir, emit: outdir)
        path("*.tsv", emit: qc_input_files)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/qc_cluster_input_files"
        cb_results_tsvs = cb_results_tsvs.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellbender__gather_qc_input: ${process_info}"
        echo "outdir: ${outdir}"
        echo "results paths: ${cb_results_tsvs}"
        045-prepare_nf_qc_cluster_input.py \\
            --cb_results_tsvs ${cb_results_tsvs}
        """
}
