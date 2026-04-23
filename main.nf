nextflow.enable.dsl=2

include { QUAST }         from './modules/qc'
include { FILTER_SAMPLES } from './modules/filter.nf'
include { AMRFINDER_SETUP } from './modules/amrfinder_setup.nf'
include { RUN_AMRFINDER } from './modules/amrfinder'
include { RUN_RGI }       from './modules/rgi'
include { RUN_ABRICATE }  from './modules/abricate'
include { RUN_RESFINDER } from './modules/resfinder'
include { HARMONIZE }     from './modules/harmonize'
include { COMPARE }       from './modules/compare'
include { RESFINDER_SETUP } from './modules/resfinder_setup.nf'
include { RGI_SETUP }       from './modules/rgi_setup.nf'

workflow {
    // 1. Inputs
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fasta)) }

    // 2. Pre-run Setup
    QUAST(samples_ch)

    QUAST.out.join(samples_ch).set { to_filter_ch }

    // 4. Run the Filter
    FILTER_SAMPLES(to_filter_ch)

    passed_samples = FILTER_SAMPLES.out.passed

    db_ch = AMRFINDER_SETUP()
    resfinder_db_ch = RESFINDER_SETUP()
    card_db_ch     = RGI_SETUP()
    
    // 3. Tool Execution with explicit output handling
    // We use 'if' inside the workflow to create the channels
    
    ch_to_mix = []
    if (params.run_amrfinder) { 
        RUN_AMRFINDER(passed_samples, db_ch.db_files)
        ch_to_mix << RUN_AMRFINDER.out
    }
    if (params.run_rgi) {
        RUN_RGI(passed_samples, card_db_ch.db_files)
        ch_to_mix << RUN_RGI.out
    }
    if (params.run_abricate) { 
        RUN_ABRICATE(passed_samples)
        ch_to_mix << RUN_ABRICATE.out
    }
    if (params.run_resfinder) {
        RUN_RESFINDER(passed_samples, resfinder_db_ch.db_files)
        ch_to_mix << RUN_RESFINDER.out
    }

    // 4. Mixing Logic
    // This takes the list of active channels and groups them
    if (ch_to_mix.size() > 0) {
        ch_to_mix[0]
            .mix( *ch_to_mix.drop(1) )
            .groupTuple()
            .set { grouped_outputs }

        HARMONIZE(grouped_outputs)
        
        // 5. Final Comparison
        COMPARE(HARMONIZE.out.collect())
    } else {
        log.info "No AMR tools selected. Skipping harmonization."
    }
}