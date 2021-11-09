#!/bin/bash

condapath=~/miniconda3
while getopts t:c:d:g: option
do
    case "${option}" in
        t) threads=${OPTARG};;
        d) download=${OPTARG};; # download address of host genome from Ensembl (e.g, http://ftp.ensembl.org/pub/release-104/fasta/callithrix_jacchus/dna/Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz)
        c) customized=${OPTARG};; # taxid for the host species
        g) gtf=${OPTARG};; # download address of host gtf from Emsenbl (e.g, http://ftp.ensembl.org/pub/release-104/gtf/callithrix_jacchus/Callithrix_jacchus.ASM275486v1.104.gtf.gz)
    esac
done

# get MTD folder place; same as Install.sh script file path (in the MTD folder)
dir=$(dirname $(readlink -f $0))
cd $dir # MTD folder place

# get conda path
condapath=$(head -n 1 $MTDIR/condaPath)
# activate MTD conda environment
source $condapath/etc/profile.d/conda.sh
conda activate MTD

# Kraken2 database building - Customized
DBNAME=kraken2DB_${customized}
mkdir -p $DBNAME
cd $DBNAME
wget -c $download #http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
unpigz *.fa.gz
mv *.fa genome_${customized}.fa
cd ..
kraken2-build --download-taxonomy --threads $threads --db $DBNAME
kraken2-build --add-to-library $DBNAME/genome_${customized}.fa --threads $threads --db $DBNAME
kraken2-build --build --threads $threads --db $DBNAME

# download host GTF
wget -c $gtf -P ref_${customized} -O ref_${customized}.gtf.gz

# Building indexes for hisat2
mkdir -p hisat2_index_${customized}
cd hisat2_index_${customized}
cp ../ref_${customized}/*.gtf.gz .
gzip -d *.gtf.gz
mv *.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
mv ../$DBNAME/genome_${customized}.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

echo "Customized host reference building is done"