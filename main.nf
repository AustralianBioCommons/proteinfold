#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/proteinfold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/proteinfold
    Website: https://nf-co.re/proteinfold
    Slack  : https://nfcore.slack.com/channels/proteinfold
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

<<<<<<< HEAD
if (params.mode == "alphafold2") {
    include { PREPARE_ALPHAFOLD2_DBS } from './subworkflows/local/prepare_alphafold2_dbs'
    include { ALPHAFOLD2             } from './workflows/alphafold2'
} else if (params.mode == "colabfold") {
    include { PREPARE_COLABFOLD_DBS } from './subworkflows/local/prepare_colabfold_dbs'
    include { COLABFOLD             } from './workflows/colabfold'
} else if (params.mode == "esmfold") {
    include { PREPARE_ESMFOLD_DBS } from './subworkflows/local/prepare_esmfold_dbs'
    include { ESMFOLD             } from './workflows/esmfold'
}

include { PIPELINE_INITIALISATION          } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { PIPELINE_COMPLETION              } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { getColabfoldAlphafold2Params     } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { getColabfoldAlphafold2ParamsPath } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLABFOLD PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

<<<<<<< HEAD
<<<<<<< HEAD
params.colabfold_alphafold2_params      = WorkflowMain.getColabfoldAlphafold2Params(params)
params.colabfold_alphafold2_params_path = WorkflowMain.getColabfoldAlphafold2ParamsPath(params)
=======
// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
>>>>>>> Template update for nf-core/tools version 2.10
=======
params.colabfold_alphafold2_params      = getColabfoldAlphafold2Params()
params.colabfold_alphafold2_params_path = getColabfoldAlphafold2ParamsPath()
>>>>>>> First iteration towards lib removal
=======
include { PROTEINFOLD  } from './workflows/proteinfold'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')
>>>>>>> Template update for nf-core/tools version 2.13

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
<<<<<<< HEAD
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_PROTEINFOLD {

    main:
    ch_multiqc  = Channel.empty()
    ch_versions = Channel.empty()

    //
    // WORKFLOW: Run alphafold2
    //
    if(params.mode == "alphafold2") {
        //
        // SUBWORKFLOW: Prepare Alphafold2 DBs
        //
        PREPARE_ALPHAFOLD2_DBS (
            params.alphafold2_db,
            params.full_dbs,
            params.bfd_path,
            params.small_bfd_path,
            params.alphafold2_params_path,
            params.mgnify_path,
            params.pdb70_path,
            params.pdb_mmcif_path,
            params.uniref30_alphafold2_path,
            params.uniref90_path,
            params.pdb_seqres_path,
            params.uniprot_path,
            params.uniprot_sprot,
            params.uniprot_trembl
        )
        ch_versions = ch_versions.mix(PREPARE_ALPHAFOLD2_DBS.out.versions)

        //
        // WORKFLOW: Run nf-core/alphafold2 workflow
        //
        ALPHAFOLD2 (
            ch_versions,
            params.full_dbs,
            params.alphafold2_mode,
            params.alphafold2_model_preset,
            PREPARE_ALPHAFOLD2_DBS.out.params,
            PREPARE_ALPHAFOLD2_DBS.out.bfd.ifEmpty([]),
            PREPARE_ALPHAFOLD2_DBS.out.small_bfd.ifEmpty([]),
            PREPARE_ALPHAFOLD2_DBS.out.mgnify,
            PREPARE_ALPHAFOLD2_DBS.out.pdb70,
            PREPARE_ALPHAFOLD2_DBS.out.pdb_mmcif,
            PREPARE_ALPHAFOLD2_DBS.out.uniref30,
            PREPARE_ALPHAFOLD2_DBS.out.uniref90,
            PREPARE_ALPHAFOLD2_DBS.out.pdb_seqres,
            PREPARE_ALPHAFOLD2_DBS.out.uniprot
        )
        ch_multiqc  = ALPHAFOLD2.out.multiqc_report
        ch_versions = ch_versions.mix(ALPHAFOLD2.out.versions)
    }

    //
    // WORKFLOW: Run colabfold
    //
    else if(params.mode == "colabfold") {
        //
        // SUBWORKFLOW: Prepare Colabfold DBs
        //
        PREPARE_COLABFOLD_DBS (
            params.colabfold_db,
            params.colabfold_server,
            params.colabfold_alphafold2_params_path,
            params.colabfold_db_path,
            params.uniref30_colabfold_path,
            params.colabfold_alphafold2_params,
            params.colabfold_db_link,
            params.uniref30_colabfold_link,
            params.create_colabfold_index
        )
        ch_versions = ch_versions.mix(PREPARE_COLABFOLD_DBS.out.versions)

        //
        // WORKFLOW: Run nf-core/colabfold workflow
        //
        COLABFOLD (
            ch_versions,
            params.colabfold_model_preset,
            PREPARE_COLABFOLD_DBS.out.params,
            PREPARE_COLABFOLD_DBS.out.colabfold_db,
            PREPARE_COLABFOLD_DBS.out.uniref30,
            params.num_recycle
        )
        ch_multiqc  = COLABFOLD.out.multiqc_report
        ch_versions = ch_versions.mix(COLABFOLD.out.versions)
    }

    //
    // WORKFLOW: Run esmfold
    //
    else if(params.mode == "esmfold") {
        //
        // SUBWORKFLOW: Prepare esmfold DBs
        //
        PREPARE_ESMFOLD_DBS (
            params.esmfold_db,
            params.esmfold_params_path,
            params.esmfold_3B_v1,
            params.esm2_t36_3B_UR50D,
            params.esm2_t36_3B_UR50D_contact_regression
        )
        ch_versions = ch_versions.mix(PREPARE_ESMFOLD_DBS.out.versions)

        //
        // WORKFLOW: Run nf-core/esmfold workflow
        //
        ESMFOLD (
            ch_versions,
            PREPARE_ESMFOLD_DBS.out.params,
            params.num_recycle
        )
        ch_multiqc  = ESMFOLD.out.multiqc_report
        ch_versions = ch_versions.mix(ESMFOLD.out.versions)
    }
    emit:
    multiqc_report = ch_multiqc  // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [version1, version2, ...]
}
=======
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PROTEINFOLD {
>>>>>>> Template update for nf-core/tools version 2.13

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    PROTEINFOLD (
        samplesheet
    )

    emit:
    multiqc_report = PROTEINFOLD.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
<<<<<<< HEAD
=======

>>>>>>> Template update for nf-core/tools version 2.13
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
<<<<<<< HEAD
        params.outdir
=======
        params.outdir,
        params.input
>>>>>>> Template update for nf-core/tools version 2.13
    )

    //
    // WORKFLOW: Run main workflow
    //
<<<<<<< HEAD
    NFCORE_PROTEINFOLD ()
=======
    NFCORE_PROTEINFOLD (
        PIPELINE_INITIALISATION.out.samplesheet
    )
>>>>>>> Template update for nf-core/tools version 2.13

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PROTEINFOLD.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
