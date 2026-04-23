process AMRFINDER_SETUP {
    conda "bioconda::ncbi-amrfinderplus"
    // Using storeDir keeps the DB across different runs
    storeDir "${params.outdir}/amrfinder_db"

    output:
    path "amrfinder_db/latest", emit: db_files

    script:
    """
    amrfinder_update --database amrfinder_db
    """
}