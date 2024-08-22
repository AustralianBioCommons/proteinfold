/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { RUN_ALPHAFOLD2      } from '../modules/local/run_alphafold2'
include { RUN_ALPHAFOLD2_MSA  } from '../modules/local/run_alphafold2_msa'
include { RUN_ALPHAFOLD2_PRED } from '../modules/local/run_alphafold2_pred'
include { EXTRACT_OUTPUTS } from '../modules/local/extract_outputs'
include { GENERATE_REPORT } from '../modules/local/generat_report'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC } from '../modules/nf-core/multiqc/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { fromSamplesheet        } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_proteinfold_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ALPHAFOLD2 {

    take:
    full_dbs                // boolean: Use full databases (otherwise reduced version)
    alphafold2_mode         //  string: Mode to run Alphafold2 in
    alphafold2_model_preset //  string: Specifies the model preset to use for Alphafold2
    ch_alphafold2_params    // channel: path(alphafold2_params)
    ch_bfd                  // channel: path(bfd)
    ch_small_bfd            // channel: path(small_bfd)
    ch_mgnify               // channel: path(mgnify)
    ch_pdb70                // channel: path(pdb70)
    ch_pdb_mmcif            // channel: path(pdb_mmcif)
    ch_uniref30             // channel: path(uniref30)
    ch_uniref90             // channel: path(uniref90)
    ch_pdb_seqres           // channel: path(pdb_seqres)
    ch_uniprot              // channel: path(uniprot)

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions = Channel.empty()

    //
    // Create input channel from input file provided through params.input
    //

    Channel
        .fromSamplesheet("input")
        .set { ch_fasta }

    if (alphafold2_model_preset != 'multimer') {
        ch_fasta
            .map {
                meta, fasta ->
                [ meta, fasta.splitFasta(file:true) ]
            }
            .transpose()
            .set { ch_fasta }
    }
    ch_alphafold_outputs = Channel.empty()

    if (alphafold2_mode == 'standard') {
        //
        // SUBWORKFLOW: Run Alphafold2 standard mode
        //
        //full_dbs.view()
        //ch_alphafold2_params.view()

        RUN_ALPHAFOLD2 (
            ch_fasta,
            full_dbs,
            alphafold2_model_preset,
            ch_alphafold2_params,
            ch_bfd,
            ch_small_bfd,
            ch_mgnify,
            ch_pdb70,
            ch_pdb_mmcif,
            ch_uniref30,
            ch_uniref90,
            ch_pdb_seqres,
            ch_uniprot
        )
        RUN_ALPHAFOLD2.out.af_out_tsv
        .map{[it[0], it[1].findAll{ it.getName().contains("_lddt_")}]}
        .set{ch_af_out_lddt}
        
        RUN_ALPHAFOLD2.out.af_out_tsv
        .map{[it[0], it[1].findAll{ it.getName().contains("_msa.tsv")}]}
        .set{ch_af_out_msa}

        RUN_ALPHAFOLD2.out.af_out_pdb.set{ch_af_out_pdb}
        
        ch_multiqc_rep = RUN_ALPHAFOLD2.out.multiqc.collect()
        ch_versions    = ch_versions.mix(RUN_ALPHAFOLD2.out.versions)
         

    } else if (alphafold2_mode == 'split_msa_prediction') {
        //
        // SUBWORKFLOW: Run Alphafold2 split mode, MSA and predicition
        //
        RUN_ALPHAFOLD2_MSA (
            ch_fasta,
            full_dbs,
            alphafold2_model_preset,
            ch_alphafold2_params,
            ch_bfd,
            ch_small_bfd,
            ch_mgnify,
            ch_pdb70,
            ch_pdb_mmcif,
            ch_uniref30,
            ch_uniref90,
            ch_pdb_seqres,
            ch_uniprot
        )
        ch_versions    = ch_versions.mix(RUN_ALPHAFOLD2_MSA.out.versions)
        
        RUN_ALPHAFOLD2_MSA.out.features
            .map{[it[0].id, it[0], it[1]]}
            .join(
         ch_fasta
            .map{[it[0].id, it[0], it[1]]}
            )
        .set{ch_af_all}

        RUN_ALPHAFOLD2_PRED (
            ch_af_all.map{[it[3], it[4]]},
            full_dbs,
            alphafold2_model_preset,
            ch_alphafold2_params,
            ch_bfd,
            ch_small_bfd,
            ch_mgnify,
            ch_pdb70,
            ch_pdb_mmcif,
            ch_uniref30,
            ch_uniref90,
            ch_pdb_seqres,
            ch_uniprot,
            ch_af_all.map{it[2]}
        )
        ch_multiqc_rep = RUN_ALPHAFOLD2_PRED.out.multiqc.map{1}.collect()
        ch_versions = ch_versions.mix(RUN_ALPHAFOLD2_PRED.out.versions)
        RUN_ALPHAFOLD2_PRED.out.af_out_tsv
        .map{[it[0], it[1].findAll{ it.getName().contains("_lddt_")}]}
        .set{ch_af_out_lddt}
        
        RUN_ALPHAFOLD2_PRED.out.af_out_tsv
        .map{[it[0], it[1].findAll{ it.getName().contains("_msa.tsv")}]}
        .set{ch_af_out_msa}

        RUN_ALPHAFOLD2_PRED.out.af_out_pdb.set{ch_af_out_pdb}
                
    }
    
    GENERATE_REPORT(
        ch_af_out_msa,
        ch_af_out_lddt,
        ch_af_out_pdb,
        Channel.fromPath("$projectDir/assets/proteinfold_template.html").first(),
        Channel.value("ALPHAFOLD2")
    )
    
  

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_proteinfold_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_report        = Channel.empty()
        ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo )   : Channel.empty()
        summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_rep)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        ch_multiqc_report = MULTIQC.out.report.toList()
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
