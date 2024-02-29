/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
<<<<<<< HEAD
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

<<<<<<< HEAD:workflows/esmfold.nf
// Validate input parameters
WorkflowEsmfold.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.esmfold_db
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input file not specified!' }
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo )   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
=======
>>>>>>> Clean esmfold
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { RUN_ESMFOLD               } from '../modules/local/run_esmfold'
include { MULTIFASTA_TO_SINGLEFASTA } from '../modules/local/multifasta_to_singlefasta'

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

workflow ESMFOLD {

    take:
    ch_versions       // channel: [ path(versions.yml) ]
    ch_esmfold_params // directory: /path/to/esmfold/params/
    ch_num_recycle    // int: Number of recycles for esmfold
    
    main:
    ch_multiqc_files = Channel.empty()

    //
    // Create input channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .set { ch_fasta }

    //
    // MODULE: Run esmfold
    //
    if (params.esmfold_model_preset != 'monomer') {
        MULTIFASTA_TO_SINGLEFASTA(
            ch_fasta
        )
        ch_versions = ch_versions.mix(MULTIFASTA_TO_SINGLEFASTA.out.versions)
        RUN_ESMFOLD(
            MULTIFASTA_TO_SINGLEFASTA.out.input_fasta,
            ch_esmfold_params,
            ch_num_recycle
        )
        ch_versions = ch_versions.mix(RUN_ESMFOLD.out.versions)
    } else {
        RUN_ESMFOLD(
            ch_fasta,
            ch_esmfold_params,
            ch_num_recycle
        )
        ch_versions = ch_versions.mix(RUN_ESMFOLD.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_proteinfold_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
<<<<<<< HEAD
    workflow_summary    = WorkflowEsmfold.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowProteinfold.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)
=======
    ch_multiqc_report                     = Channel.empty()
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo )   : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_methods_description))
>>>>>>> Clean esmfold

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(RUN_ESMFOLD.out.multiqc.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_multiqc_report = MULTIQC.out.report.toList()
    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
