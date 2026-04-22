process RUN_AMRFINDER {
    tag "$sample"
    publishDir "${params.outdir}/amrfinder/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_amrfinder.tsv")

    script:
    """
    amrfinder \
        --nucleotide ${fasta} \
        --output ${sample}_amrfinder.tsv \
        --ident-min ${params.amrfinder_ident ?: 90} \
        --coverage-min ${params.amrfinder_cov ?: 50} \
        --threads 4
    """
}