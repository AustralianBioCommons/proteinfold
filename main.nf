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

params.colabfold_alphafold2_params      = WorkflowMain.getColabfoldAlphafold2Params(params)
params.colabfold_alphafold2_params_path = WorkflowMain.getColabfoldAlphafold2ParamsPath(params)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_PROTEINFOLD {

    main:
    ch_multiqc  = Channel.empty()
    ch_versions = Channel.empty()

    //
    // WORKFLOW: Run alphafold2
    //
    if(params.mode.toLowerCase().split(",").contains("alphafold2")) {
        //
        // SUBWORKFLOW: Prepare Alphafold2 DBs
        //
        // WORKFLOW: Run nf-core/alphafold2 workflow
        //

        ch_params         = Channel.fromPath( params.alphafold2_params_path )
        ch_mgnify         = Channel.fromPath( params.mgnify_path )
        ch_pdb70          = Channel.fromPath( params.pdb70_path, type: 'dir' )
        ch_mmcif_files    = Channel.fromPath( params.pdb_mmcif_path, type: 'dir' )
        ch_mmcif_obsolete = Channel.fromPath( params.pdb_mmcif_path, type: 'file' )
        ch_mmcif          = ch_mmcif_files.mix(ch_mmcif_obsolete)
        ch_uniref30       = Channel.fromPath( params.uniref30_alphafold2_path, type: 'any' )
        ch_uniref90       = Channel.fromPath( params.uniref90_path )
        ch_pdb_seqres     = Channel.fromPath( params.pdb_seqres_path )
        ch_uniprot        = Channel.fromPath( params.uniprot_path )
        ch_small_bfd = Channel.fromPath( params.small_bfd_path)
        ch_bfd = Channel.fromPath( params.bfd_path)
        
        
        
        
        ALPHAFOLD2 (
            params.full_dbs,
            params.alphafold2_mode,
            params.alphafold2_model_preset,
            ch_params.toList(),
            ch_bfd.ifEmpty([]).first(),
            ch_small_bfd.ifEmpty([]).first(),
            ch_mgnify.first(),
            ch_pdb70.first(),
            ch_mmcif.toList(),
            ch_uniref30.toList(),
            ch_uniref90.first(),
            ch_pdb_seqres.first(),
            ch_uniprot.first()
        )
        ch_multiqc  = ALPHAFOLD2.out.multiqc_report
        ch_versions = ch_versions.mix(ALPHAFOLD2.out.versions)
    }

    //
    // WORKFLOW: Run colabfold
    //
    if(params.mode.toLowerCase().split(",").contains("colabfold")) {
        //
        // SUBWORKFLOW: Prepare Colabfold DBs
        //
        PREPARE_COLABFOLD_DBS (
            params.colabfold_db,
            params.colabfold_server,
            params.colabfold_alphafold2_params_path,
            params.colabfold_db_path,
            params.uniref30_colabfold_path,
            params.colabfold_alphafold2_params_link,
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
            PREPARE_COLABFOLD_DBS.out.params.first(),
            PREPARE_COLABFOLD_DBS.out.colabfold_db.first(),
            PREPARE_COLABFOLD_DBS.out.uniref30.first(),
            params.num_recycles_colabfold
        )
        ch_multiqc  = COLABFOLD.out.multiqc_report
        ch_versions = ch_versions.mix(COLABFOLD.out.versions)
    }

    //
    // WORKFLOW: Run esmfold
    //
    if(params.mode.toLowerCase().split(",").contains("esmfold")) {
        //
        // SUBWORKFLOW: Prepare esmfold DBs
        //
        /*PREPARE_ESMFOLD_DBS (
            params.esmfold_db,
            params.esmfold_params_path,
            params.esmfold_3B_v1,
            params.esm2_t36_3B_UR50D,
            params.esm2_t36_3B_UR50D_contact_regression
        )*/
        
        //ch_versions = ch_versions.mix(PREPARE_ESMFOLD_DBS.out.versions)

        //
        // WORKFLOW: Run nf-core/esmfold workflow
        //
        ESMFOLD (
            ch_versions,
            Channel.fromPath(params.esmfold_params_path),
            params.num_recycles_esmfold
        )
        ch_multiqc  = ESMFOLD.out.multiqc_report
        ch_versions = ch_versions.mix(ESMFOLD.out.versions)
    }
    emit:
    multiqc_report = ch_multiqc  // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [version1, version2, ...]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PROTEINFOLD ()

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
