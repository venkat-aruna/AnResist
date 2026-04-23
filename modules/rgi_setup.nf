process RGI_SETUP {
    conda "bioconda::rgi=6.0.3"
    storeDir "${params.outdir}/card_db"

    output:
    path "card.json", emit: db_files

    script:
    """
    wget https://card.mcmaster.ca/latest/data -O data.tar.bz2
    tar -xvf data.tar.bz2
    """
}