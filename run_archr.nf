#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


// params.rpath = "/software/isg/languages/R/4.4.0/exec/bin/Rscript"
params.sample_list_path = "./data/test_snapatac2_files.csv"

process archr_per_sample {
    label "archr"
    publishDir "arrowfiles", mode: 'copy'

    input:
        val samplename
    
    output:
        path "${samplename}_genescore.mtx"
        path "${samplename}_features.tsv"
        path "${samplename}_barcodes.tsv"

    script:
    """
    Rscript ${workflow.projectDir}/scripts/archr_process.R ${workflow.projectDir}/data/fragment_files/${samplename}/fragments.tsv.gz ${samplename}
    """
}

workflow {
    sample_ch = Channel.fromPath(params.sample_list_path)
                        .splitText()
                        .map { it.trim() }
                        .filter { it }
                        .view { "Process ${it}" }
                        | archr_per_sample

}