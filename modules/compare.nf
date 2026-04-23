process COMPARE {
    conda "$projectDir/environment.yml"
    publishDir "${params.outdir}/comparison", mode: 'copy'

    input:
    path(unified_tsvs)

    output:
    path("comparison/")

    script:
    """
    python ${projectDir}/src/compare.py \
        --input ${unified_tsvs} \
        --output_dir comparison/
    """
}