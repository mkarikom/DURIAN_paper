#!/bin/bash
module load salmon

transfolder=/dfs5/bio/mkarikom/code/DTMwork/slurm/He/human_transcriptome
mkdir -p $transfolder

# get the data
wget -P $transfolder \
ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/\
Homo_sapiens.GRCh38.cdna.all.fa.gz

# Decompress the FASTA file
gunzip -c $transfolder/Homo_sapiens.GRCh38.cdna.all.fa.gz > \
$transfolder/Homo_sapiens.GRCh38.cdna.all.fa

indexfolder=/dfs5/bio/mkarikom/code/DTMwork/slurm/He/salmon_index
mkdir -p $indexfolder

salmon index \
-t $transfolder/Homo_sapiens.GRCh38.cdna.all.fa \
-i $indexfolder \
-k 31