#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.atac_names_file = "./data/test_snapatac2_files.csv"

process test {
    publishDir "data/test_results", mode: 'copy'

    input:
    val name

    output:
    path "${name}.csv"

    script:
    """
    python ${workflow.projectDir}/scripts/test.py ${workflow.projectDir}/data/fragment_files/${name}
    echo ${name} > ${name}.csv
    """
}


workflow {

    atac_names = Channel
        .fromPath(params.atac_names_file)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .view { "Processing line: $it" } 
        | test

}

