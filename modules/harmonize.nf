// modules/harmonize.nf
process HARMONIZE {
    tag "$sample"
    publishDir "${params.outdir}/harmonized", mode: 'copy'

    input:
    tuple val(sample), path(amrfinder), path(rgi), path(abricate),
          path(resfinder), path(staramr), path(fargene)

    output:
    path("${sample}_unified.tsv")

    script:
    """
    python ${projectDir}/src/parser.py \
        --sample ${sample} \
        --amrfinder ${amrfinder} \
        --rgi ${rgi} \
        --abricate ${abricate} \
        --resfinder ${resfinder} \
        --staramr ${staramr} \
        --fargene ${fargene} \
        --output ${sample}_unified.tsv
    """
}