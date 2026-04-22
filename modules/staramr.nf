process RUN_STARAMR {
    tag "$sample"
    conda "bioconda::staramr"
    publishDir "${params.outdir}/staramr/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_staramr.tsv")

    script:
    """
    staramr search \
        --output-dir staramr_out \
        ${fasta}

    cp staramr_out/resfinder.tsv ${sample}_staramr.tsv
    """
}