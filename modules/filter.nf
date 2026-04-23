process FILTER_SAMPLES {
    tag "$sample"
    conda "conda-forge::python=3.10"

    input:
    tuple val(sample), path(quast_dir), path(fasta)

    output:
    // This channel only emits samples that pass
    tuple val(sample), path(fasta), emit: passed
    // Optional: emit failed samples for a log
    tuple val(sample), path("fail_reason.txt"), emit: failed, optional: true

    script:
    """
    python <<EOF
    import csv
    import os

    # Load parameters from Nextflow environment
    min_size = int('${params.min_genome_size}')
    max_contigs = int('${params.min_contigs}') # Assuming this is a max limit

    passed = False
    fail_msg = ""

    with open('${quast_dir}/transposed_report.tsv', 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        row = next(reader)
        genome_size = int(row['Total length'])
        contigs = int(row['# contigs'])

        if genome_size >= min_size:
            passed = True
        else:
            fail_msg = f"Size {genome_size} < {min_size}"

    if passed:
        # Create a symbolic link or just let Nextflow handle the 'passed' output
        os.system('ln -s ${fasta} passed_fasta.fasta')
    else:
        with open('fail_reason.txt', 'w') as f:
            f.write(fail_msg)
    EOF
    """
}