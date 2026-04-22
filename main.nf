nextflow.enable.dsl=2

include { QUAST }         from './modules/qc'
include { RUN_AMRFINDER } from './modules/amrfinder'
include { RUN_RGI }       from './modules/rgi'
include { RUN_ABRICATE }  from './modules/abricate'
include { RUN_RESFINDER } from './modules/resfinder'
include { RUN_STARAMR }   from './modules/staramr'
include { RUN_FARGENE }   from './modules/fargene'
include { HARMONIZE }     from './modules/harmonize'
include { COMPARE }       from './modules/compare'

workflow {

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fasta)) }
        .set { samples }

    QUAST(samples)

    // run all tools unconditionally — use params to filter inputs instead
    if (params.run_amrfinder) { RUN_AMRFINDER(samples) }
    if (params.run_rgi)       { RUN_RGI(samples) }
    if (params.run_abricate)  { RUN_ABRICATE(samples) }
    if (params.run_resfinder) { RUN_RESFINDER(samples) }
    if (params.run_staramr)   { RUN_STARAMR(samples) }
    if (params.run_fargene)   { RUN_FARGENE(samples) }

    // collect outputs — use ifEmpty to handle disabled tools gracefully
    amrfinder_out = params.run_amrfinder ? RUN_AMRFINDER.out : Channel.empty()
    rgi_out       = params.run_rgi       ? RUN_RGI.out       : Channel.empty()
    abricate_out  = params.run_abricate  ? RUN_ABRICATE.out  : Channel.empty()
    resfinder_out = params.run_resfinder ? RUN_RESFINDER.out : Channel.empty()
    staramr_out   = params.run_staramr   ? RUN_STARAMR.out   : Channel.empty()
    fargene_out   = params.run_fargene   ? RUN_FARGENE.out   : Channel.empty()

    // mix all outputs and group by sample name
    amrfinder_out
        .mix(rgi_out, abricate_out, resfinder_out, staramr_out, fargene_out)
        .groupTuple()
        .set { all_tool_outputs }

    HARMONIZE(all_tool_outputs)
    COMPARE(HARMONIZE.out.collect())
}