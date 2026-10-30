//
// Subworkflow with functionality specific to the duplex-seq-pipeline pipeline
//

include { validateParameters ; paramsSummaryLog ; samplesheetToList } from 'plugin/nf-schema'

workflow PIPELINE_INITIALISATION {
    take:
    validate_params // boolean: validate parameters against nextflow_schema.json
    input // string: path to the samplesheet

    main:

    if (validate_params) {
        validateParameters()
    }

    log.info(paramsSummaryLog(workflow))

    // [ meta, fastq_1, fastq_2, umi_file ] per samplesheet row. Relative paths in the
    // samplesheet are resolved against the samplesheet's directory.
    def samplesheet_dir = file(input).parent
    ch_samplesheet = channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2, umi_file ->
            def umis = umi_file ? resolveSamplesheetPath(umi_file, samplesheet_dir) : (params.umis ? file(params.umis, checkIfExists: true) : [])
            [
                meta + [single_end: false, read_structure: meta.read_structure ?: params.read_structure],
                [resolveSamplesheetPath(fastq_1, samplesheet_dir), resolveSamplesheetPath(fastq_2, samplesheet_dir)],
                umis,
            ]
        }

    emit:
    samplesheet = ch_samplesheet // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ], path(umi_file) or [] ]
}

//
// Resolve a samplesheet path: absolute paths and URLs are used as they are
//
def resolveSamplesheetPath(path, samplesheet_dir) {
    if (path.startsWith('/') || path.contains('://')) {
        return file(path, checkIfExists: true)
    }
    return file(samplesheet_dir.resolve(path), checkIfExists: true)
}

//
// Convert the `versions` topic into a YAML file of tool versions
//
def softwareVersionsToYAML(ch_topic_versions) {
    def workflow_versions = "\"Workflow\":\n    ${workflow.manifest.name}: ${workflow.manifest.version}\n    Nextflow: ${workflow.nextflow.version}\n".toString()
    return ch_topic_versions
        .unique()
        .map { proc, tool, version -> [proc.tokenize(':')[-1], "${tool}: ${version.trim()}".toString()] }
        .groupTuple()
        .map { proc, tools -> ("\"${proc}\":\n" + tools.sort().collect { tool -> "    ${tool}" }.join('\n') + '\n').toString() }
        .mix(channel.value(workflow_versions))
        .collectFile(name: 'software_versions.yml', sort: true, newLine: false, storeDir: "${params.outdir}/pipeline_info")
}
