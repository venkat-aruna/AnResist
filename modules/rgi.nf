process RUN_RGI {
    tag "$sample"
    publishDir "${params.outdir}/rgi/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_rgi.txt")

    script:
    def loose = params.rgi_loose ? "--include_loose" : ""
    """
    rgi main \
        --input_sequence ${fasta} \
        --output_file ${sample}_rgi \
        --input_type contig \
        --num_threads 4 \
        --clean \
        ${loose}

    mv ${sample}_rgi.txt ${sample}_rgi.txt
    """
}