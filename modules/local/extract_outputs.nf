process EXTRACT_OUTPUTS {
    tag "$id"
    label 'process_single'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ALPHAFOLD2 module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    container "nf-core/proteinfold_alphafold2_standard:dev"

    input:
    tuple val(id), path(pkl_files)
    
    output:
    tuple val(id), path ("*_msa.tsv"), optional: true, emit: msa_info
    tuple val(id), path ("*_lddt_*.tsv"), optional: true, emit: lddt_info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    //#for pkl_file in [\"${"\", \"".join(pkl_files)}\"]:
    """
    #!/usr/bin/env python
    import pickle
    import os, sys
    for pkl_file in [\"${pkl_files.join("\", \"")}\"]:
        dict_data = pickle.load(open(pkl_file,'rb'))
        if pkl_file.endswith("features.pkl"):
            with open ("${id}_msa.tsv", "w") as out_f:
                for val in dict_data['msa']:
                    out_f.write("\\t".join([str(x) for x in val]) + "\\n")
        else:
            model_id = os.path.basename(pkl_file).replace("result_model_", "").replace("_pred_0.pkl", "") 
            with open (f"${id}_lddt_{model_id}.tsv", "w") as out_f:
                out_f.write("\\t".join([str(x) for x in dict_data['plddt']]) + "\\n")
    with open ("versions.yml", "w") as version_file:
        version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))
    """
}
