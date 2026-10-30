#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    duplex-seq-pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/mcmero/duplex-seq-pipeline
----------------------------------------------------------------------------------------
*/

include { DUPLEXSEQ } from './workflows/duplexseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_pipeline/main'

workflow {

    PIPELINE_INITIALISATION(
        params.validate_params,
        params.input,
    )

    DUPLEXSEQ(PIPELINE_INITIALISATION.out.samplesheet)
}
