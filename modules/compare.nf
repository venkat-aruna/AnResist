// modules/compare.nf
process COMPARE {
    conda "$projectDir/environment.yml"
    publishDir "${params.outdir}/comparison", mode: 'copy'
    
    input:
    path(unified_files)

    output:
    path("comparison_results/")

    script:
    """
    python ${projectDir}/src/compare.py \
        --input ${unified_files} \
        --output comparison_results/
    """
}