process GENERATE_REPORT {
    tag "${meta.id}"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.21--pyhdfd78af_0' :
        'biocontainers/multiqc:1.21--pyhdfd78af_0' }"

    input:
    tuple val(meta_msa), path(msa)
    tuple val(meta), path(lddt)
    tuple val(meta), path(pdb)
    path(template)
    val(output_type)
    
    output:
    tuple val(meta), path ("*.html"), emit: report
    tuple val(meta), path ("*.png"), emit: images
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    
    """
    generat_plots_2.py --type ${output_type} --msa ${msa} --plddt ${lddt.join(' ')} --pdb ${pdb.join(' ')} --html_template ${template} --output_dir ./ --name ${meta.id}
    """
}
