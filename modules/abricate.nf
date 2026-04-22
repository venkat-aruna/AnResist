process RUN_ABRICATE {
    tag "$sample"
    publishDir "${params.outdir}/abricate/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_abricate.tsv")

    script:
    """
    abricate \
        --db ${params.abricate_db ?: 'ncbi'} \
        --minid ${params.abricate_minid ?: 80} \
        --mincov ${params.abricate_mincov ?: 80} \
        ${fasta} > ${sample}_abricate.tsv
    """
}