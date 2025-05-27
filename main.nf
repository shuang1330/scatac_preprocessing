#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters must be declared at the top
params.concated_path     = "concated.h5ads"
params.peakmat_path      = "peakmat.h5ad"
params.atac_names_file   = "./data/snapatac2_files.csv"
params.fragments_dir = "./data/fragment_files"
params.snapatac_dir = "./data/snapatac2_files"

process process_per_fragment {
    label "preprocess"
    publishDir "data/snapatac2_files", mode: 'copy'

    input:
    val name

    output:
    path "${name}.h5ad", emit: fragment_files

    script:
    """
    python ${workflow.projectDir}/scripts/snapatac2_per_fragment.py \
        --fragments_dir ${workflow.projectDir}/${params.fragments_dir} \
        --name ${name}
    """
}

process concat_peakcalling {
    label "peakcalling"

    publishDir "data", mode: 'copy'

    input:
    val snapatac_dir
    val concated_path
    val peakmat_path

    output: 
    path "${concated_path}", emit: anndataset
    path "${peakmat_path}", emit: peakmat

    script:
    """
    python ${workflow.projectDir}/scripts/snapatac2_peak_calling.py \
        --snapatac_dir ${workflow.projectDir}/${snapatac_dir} \
        --concated_filename ${concated_path} \
        --peakmat_savepath ${peakmat_path}
    """
}

workflow {

    atac_names = Channel
        .fromPath(params.atac_names_file)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .view { "Processing line: $it" } 
        | process_per_fragment

    process_per_fragment.out.fragment_files
        .collect()
        .map { _ -> params.snapatac_dir }
        .set { snapatac_dir_ready }

    concat_peakcalling(
        snapatac_dir_ready,
        params.concated_path,
        params.peakmat_path
    )
}
