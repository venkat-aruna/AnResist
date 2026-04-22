nextflow.enable.dsl=2

// import modules
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

    // read samplesheet → emit (sample_name, fasta_path) tuples
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fasta)) }
        .set { samples }

    // QC always runs
    QUAST(samples)

    // AMR tools — each gated by its toggle in params
    amrfinder_out = params.run_amrfinder ? RUN_AMRFINDER(samples) : Channel.empty()
    rgi_out       = params.run_rgi       ? RUN_RGI(samples)       : Channel.empty()
    abricate_out  = params.run_abricate  ? RUN_ABRICATE(samples)  : Channel.empty()
    resfinder_out = params.run_resfinder ? RUN_RESFINDER(samples) : Channel.empty()
    staramr_out   = params.run_staramr   ? RUN_STARAMR(samples)   : Channel.empty()
    fargene_out   = params.run_fargene   ? RUN_FARGENE(samples)   : Channel.empty()

    // collect all tool outputs per sample, then harmonize
    amrfinder_out
        .join(rgi_out)
        .join(abricate_out)
        .join(resfinder_out)
        .join(staramr_out)
        .join(fargene_out)
        .set { all_tool_outputs }

    HARMONIZE(all_tool_outputs)

    // collect all harmonized outputs, then compare
    COMPARE(HARMONIZE.out.collect())
}