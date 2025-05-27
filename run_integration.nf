#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.model_name = 'poissonvi'
params.data_name = 'poissonvi.h5ad'
params.params_name = 'integration_params.xml'

process run_integration {
    label "integration"

    publishDir "data/integration/${model_name}/", mode: 'copy'

    input:
    val model_name
    val data_name
    val params_name

    output: 
    path "model.pt", emit: model
    path "adat.h5ad", emit: adat

    script:
    """
    python ${workflow.projectDir}/scripts/run_intergration.py \
        --data_path ${workflow.projectDir}/data/${data_name} \
        --model_name ${model_name} \
        --params_path ${workflow.projectDir}/params/${params_name}
    """
}


workflow {
    
    run_integration(
        params.model_name,
        params.data_name,
        params.params_name
    )

}
