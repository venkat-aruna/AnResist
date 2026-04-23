process AMRFINDER_SETUP {
    conda "bioconda::ncbi-amrfinderplus"
    storeDir "${params.outdir}/amrfinder_db"

    output:
    path "amrfinder_db/2026-03-24.1", emit: db_files

    script:
    """
    amrfinder_update --database amrfinder_db
    """
}