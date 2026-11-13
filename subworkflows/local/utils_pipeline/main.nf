//
// Subworkflow with functionality specific to the duplex-seq-pipeline pipeline
//

include { validateParameters ; paramsSummaryLog ; samplesheetToList } from 'plugin/nf-schema'

workflow PIPELINE_INITIALISATION {
    take:
    validate_params // boolean: validate parameters against nextflow_schema.json
    input // string: path to the samplesheet
    merge_lanes // boolean: the rows of a sample are its lanes and are merged

    main:

    if (validate_params) {
        validateParameters()
    }

    log.info(paramsSummaryLog(workflow))

    // One entry per sample: [ meta, fastq files, umi_file ]. Relative paths in the samplesheet are
    // resolved against the samplesheet's directory. The rows are checked here, before anything runs.
    def samplesheet_dir = file(input).parent
    def samples = groupLanes(samplesheetToList(input, "${projectDir}/assets/schema_input.json"), merge_lanes, samplesheet_dir)
    ch_samplesheet = channel.fromList(samples)

    emit:
    samplesheet = ch_samplesheet
}

//
// The lane and the read of a FASTQ file name such as sample_L001_R1.fastq.gz (the bcl2fastq form
// sample_S1_L001_R1_001.fastq.gz works too): [ lane: 1, read: 1 ], or null if the name has no lane
//
def laneOfFastq(fastq) {
    def matcher = (fastq.tokenize('/')[-1] =~ /_L(\d+)_R([12])(?:_\d+)?\.f(?:ast)?q\.gz$/)
    return matcher.find() ? [lane: matcher.group(1).toInteger(), read: matcher.group(2).toInteger()] : null
}

//
// A lane number as in the file names: 1 is 001
//
def laneName(lane) {
    return String.format('%03d', lane)
}

//
// Check the rows of a sample that is on several lines of the samplesheet (its lanes) and return them
// in lane order
//
def sortLanes(sample, sample_rows) {
    def lanes = sample_rows.collect { row ->
        def fastq_1 = row[1]
        def fastq_2 = row[2]
        def lane_1 = laneOfFastq(fastq_1)
        def lane_2 = laneOfFastq(fastq_2)
        if (!lane_1 || lane_1.read != 1) {
            error("Sample '${sample}': cannot read the lane of read 1 from '${fastq_1}'. With --merge_lanes the names of the files must end like _L001_R1.fastq.gz")
        }
        if (!lane_2 || lane_2.read != 2) {
            error("Sample '${sample}': cannot read the lane of read 2 from '${fastq_2}'. With --merge_lanes the names of the files must end like _L001_R2.fastq.gz")
        }
        if (lane_1.lane != lane_2.lane) {
            error("Sample '${sample}': read 1 and read 2 of one row are from different lanes (L${laneName(lane_1.lane)} in '${fastq_1}', L${laneName(lane_2.lane)} in '${fastq_2}')")
        }
        lane_1.lane
    }

    def repeated = lanes.countBy { lane -> lane }.findAll { _lane, count -> count > 1 }.keySet()
    if (repeated) {
        error("Sample '${sample}': lane ${repeated.collect { lane -> 'L' + laneName(lane) }.join(', ')} is on more than one row of the samplesheet")
    }
    if (sample_rows.collect { row -> row[0].read_structure }.unique().size() > 1) {
        error("Sample '${sample}': the rows of its lanes have different read structures")
    }
    if (sample_rows.collect { row -> row[3] ?: '' }.unique().size() > 1) {
        error("Sample '${sample}': the rows of its lanes have different UMI files")
    }

    // numeric order, so that L2 comes before L10
    return [sample_rows, lanes].transpose().sort { pair -> pair[1] }.collect { pair -> pair[0] }
}

//
// Turn the rows of the samplesheet ([ meta, fastq_1, fastq_2, umi_file ]) into one entry per sample.
// With merge_lanes the rows of a sample are its lanes: they are put in lane order and their files are
// listed one after the other (R1, R2, R1, R2, ...), as CAT_FASTQ expects them, and meta.n_lanes is the
// number of lanes. Without merge_lanes a sample on several rows is an error.
//
def groupLanes(rows, merge_lanes, samplesheet_dir) {
    def samples = rows.groupBy { row -> row[0].id }
    return samples.collect { sample, sample_rows ->
        if (sample_rows.size() > 1 && !merge_lanes) {
            error("Sample '${sample}' is on ${sample_rows.size()} rows of the samplesheet. If these are the lanes of one sample, use --merge_lanes; otherwise give every sample a unique name.")
        }
        def lane_rows = sample_rows.size() > 1 ? sortLanes(sample, sample_rows) : sample_rows
        def first_row = lane_rows[0]
        def umi_file = first_row[3]
        def umis = umi_file ? resolveSamplesheetPath(umi_file, samplesheet_dir) : (params.umis ? file(params.umis, checkIfExists: true) : [])
        def reads = lane_rows.collectMany { row -> [resolveSamplesheetPath(row[1], samplesheet_dir), resolveSamplesheetPath(row[2], samplesheet_dir)] }
        [
            first_row[0] + [single_end: false, read_structure: first_row[0].read_structure ?: params.read_structure, n_lanes: lane_rows.size()],
            reads,
            umis,
        ]
    }
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
