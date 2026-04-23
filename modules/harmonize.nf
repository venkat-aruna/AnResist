process HARMONIZE {
    tag "$sample"
    conda "$projectDir/environment.yml"
    publishDir "${params.outdir}/harmonized", mode: 'copy'

    input:
    tuple val(sample), path(tool_files) // tool_files will be a list of TSVs

    output:
    path("${sample}_unified.tsv")

    script:
    // We pass the list of files directly to the script
    // Nextflow expands ${tool_files} into a space-separated string of filenames
    """
    python ${projectDir}/src/harmonize.py \
        --sample ${sample} \
        --input_files ${tool_files} \
        --output ${sample}_unified.tsv
    """
}