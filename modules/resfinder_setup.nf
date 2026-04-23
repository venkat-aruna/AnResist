process RESFINDER_SETUP {
    conda "bioconda::resfinder"
    storeDir "${params.outdir}/resfinder_db"

    output:
    path "resfinder_db", emit: db_files

    script:
    """
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db
    """
}