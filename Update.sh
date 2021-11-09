#!/bin/bash

condapath=~/miniconda3
while getopts t: option
do
    case "${option}" in
        t) threads=${OPTARG};;
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

# update viruses4kraken.fa
rm virushostdb.genomic.fna viruses4kraken.fa
wget -c https://www.genome.jp/ftp/db/virushostdb/virushostdb.genomic.fna.gz
unpigz virushostdb.genomic.fna.gz
cat Installation/M33262_SIVMM239.fa virushostdb.genomic.fna > viruses4kraken.fa

# Kraken2 database update/re-building - Microbiome
rm -rf kraken2DB_micro
DBNAME=kraken2DB_micro
kraken2-build --download-taxonomy --threads $threads --db $DBNAME
kraken2-build --download-library archaea --threads $threads --db $DBNAME
kraken2-build --download-library bacteria --threads $threads --db $DBNAME
kraken2-build --download-library protozoa --threads $threads --db $DBNAME
kraken2-build --download-library fungi --threads $threads --db $DBNAME
kraken2-build --download-library plasmid --threads $threads --db $DBNAME
kraken2-build --download-library UniVec_Core --threads $threads --db $DBNAME
kraken2-build --add-to-library viruses4kraken.fa --threads $threads --db $DBNAME
kraken2-build --build --threads $threads --db $DBNAME

# Bracken database update/re-building
bracken-build -d $dir/kraken2DB_micro -t $threads -l 75

echo "update done"