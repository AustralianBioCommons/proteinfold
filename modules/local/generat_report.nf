process GENERATE_REPORT {
    tag "$id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.21--pyhdfd78af_0' :
        'biocontainers/multiqc:1.21--pyhdfd78af_0' }"

    input:
    tuple val(id), path(msa)
    tuple val(id), path(lddt)
    tuple val(id), path(pdb)
    path(template)
    path(script)
    output:
    tuple val(id), path ("*.html"), emit: report
    tuple val(id), path ("*.png"), emit: images
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    
    """
    #export MPLCONFIGDIR=\$PBS_JOBFS
    python ./generat_plots_2.py --msa ${msa} --plddt ${lddt.join(' ')} --pdb ${pdb.join(' ')} --html_template ${template} --output_dir ./ --name ${id} || true
    """
}
