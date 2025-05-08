#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


params.rpath = "/software/isg/languages/R/4.4.0/exec/bin/Rscript"
params.rscriptpath = "/lustre/scratch126/cellgen/team205/abf/new_atac_samples/archr_analysis/create_arrowfile.r"
params.sample_list_path = "fragment_filelist.csv"
params.savedir = "/lustre/scratch126/cellgen/team205/abf/new_atac_samples/arrowfiles"

process archr_per_sample {

    publishDir "${params.savedir}", mode: 'copy'

    input:
        path fragment_path
        val samplename
    
    output:
        path "${samplename}/${samplename}_genescore.mtx"
        path "${samplename}/${samplename}_features.tsv"
        path "${samplename}/${samplename}_barcodes.tsv"

    script:
    """
    ${params.rpath} ${params.rscriptpath} ${fragment_path} ${samplename}
    """
}

workflow {
    sample_ch = Channel.fromPath(params.sample_list_path)
                        .splitCSv(sep: ',', header: false)
                        .map { row -> [row[0].strim(), row[1].trim()] }
                        .filter { row[0], row[1] }
                        .view { "Process ${row[1]}" }
                        | archr_per_sample

}