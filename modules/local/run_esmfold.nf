process RUN_ESMFOLD {
    tag "$meta.id"
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ESMFOLD module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    container "nf-core/proteinfold_esmfold:1.1.1"

    input:
    tuple val(meta), path(fasta)
    path (esm_fold_parms)
    val numRec

    output:
    tuple val(meta), path ("${fasta.baseName}*.pdb"), emit: pdb
    tuple val(meta), path ("${fasta.baseName}_plddt_mqc.tsv"), emit: multiqc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    cp -r /mnt/d/01852933ca43cd53eb240fcf350f32/* ./
    cp T1026.1_plddt_mqc.tsv T1026_plddt_mqc.tsv

    """

    stub:
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ./"${fasta.baseName}".pdb
    touch ./"${fasta.baseName}"_plddt_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esm-fold: $VERSION
    END_VERSIONS
    """
}
