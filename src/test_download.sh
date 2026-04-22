#!/usr/bin/bash

cd data/test
mkdir -p fastas

# Download all 10 directly via wget from NCBI FTP
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz" -O fastas/M_tuberculosis.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz" -O fastas/B_subtilis.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/505/GCF_000011505.1_ASM1150v1/GCF_000011505.1_ASM1150v1_genomic.fna.gz" -O fastas/S_aureus_MRSA.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz" -O fastas/K_pneumoniae.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz" -O fastas/P_aeruginosa.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/085/GCF_000012085.1_ASM1208v1/GCF_000012085.1_ASM1208v1_genomic.fna.gz" -O fastas/A_baumannii.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/172/575/GCF_000172575.2_ASM17257v2/GCF_000172575.2_ASM17257v2_genomic.fna.gz" -O fastas/E_faecium.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz" -O fastas/S_enterica.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/845/GCF_000006845.1_ASM684v1/GCF_000006845.1_ASM684v1_genomic.fna.gz" -O fastas/N_gonorrhoeae.fasta.gz
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/205/GCF_000009205.2_ASM920v2/GCF_000009205.2_ASM920v2_genomic.fna.gz" -O fastas/C_difficile.fasta.gz

# Decompress all
gunzip fastas/*.gz

# Verify
ls -lh fastas/