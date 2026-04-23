process RUN_RGI {
    tag "$sample"
    conda "bioconda::rgi=6.0.3"
    publishDir "${params.outdir}/rgi/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)
    path card_json

    output:
    tuple val(sample), path("${sample}_rgi.txt")

    script:
    def loose = params.rgi_loose ? "--include_loose" : ""
    """
    rgi load -i ${card_json} --local
    rgi main \
        --input_sequence ${fasta} \
        --output_file ${sample}_rgi \
        --input_type contig \
        --num_threads 4 \
        --clean \
        --local \
        ${loose}
    """
}