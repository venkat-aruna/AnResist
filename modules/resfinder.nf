process RUN_RESFINDER {
    tag "$sample"
    conda "bioconda::resfinder"
    publishDir "${params.outdir}/resfinder/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_resfinder.tsv")

    script:
    """
    run_resfinder.py \
        --inputfasta ${fasta} \
        --outputPath resfinder_out \
        --acquired \
        --threshold ${params.resfinder_threshold ?: 0.9} \
        --min_cov ${params.resfinder_cov ?: 0.6}

    cp resfinder_out/ResFinder_results_tab.txt ${sample}_resfinder.tsv
    """
}