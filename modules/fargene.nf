process RUN_FARGENE {
    tag "$sample"
    conda "bioconda::fargene"
    publishDir "${params.outdir}/fargene/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_fargene.tsv")

    script:
    """
    fargene \
        --infiles ${fasta} \
        --hmm-model ${params.fargene_hmm ?: 'class_A'} \
        --output fargene_out \
        --threads 4

    # fargene outputs a summary — convert to tsv
    cp fargene_out/results_summary.txt ${sample}_fargene.tsv
    """
}