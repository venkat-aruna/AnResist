process QUAST {
    tag "$sample"
    publishDir "${params.outdir}/qc/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("quast_results/")

    script:
    """
    quast.py ${fasta} \
        --output-dir quast_results \
        --min-contig 500 \
        --threads 4
    """
}