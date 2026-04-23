process RUN_AMRFINDER {
    tag "$sample"
    conda "bioconda::ncbi-amrfinderplus"
    publishDir "${params.outdir}/amrfinder/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)
    path amrfinder_db  // This is the new input dependency

    output:
    tuple val(sample), path("${sample}_amrfinder.tsv")

    script:
    """
    amrfinder \
        --nucleotide ${fasta} \
        --database ${amrfinder_db} \
        --output ${sample}_amrfinder.tsv \
        --ident_min ${params.amrfinder_ident} \
        --coverage_min ${params.amrfinder_cov} \
        --threads 4
    """
}