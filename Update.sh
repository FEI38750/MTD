#!/bin/bash

kmer="" # --kmer-len in kraken2-build
min_l="" # --minimizer-len in kraken2-build
min_s="" # --minimizer-spaces in kraken2-build
read_len=75 # the read length in bracken-build

condapath=~/miniconda3
while getopts t:k:m:s:r: option
do
    case "${option}" in
        t) threads=${OPTARG};;
        k) kmer=${OPTARG};; # --kmer-len in kraken2-build
        m) min_l=${OPTARG};; # --minimizer-len in kraken2-build
        s) min_s=${OPTARG};; # --minimizer-spaces in kraken2-build
        r) read_len=${OPTARG};; # the read length in bracken-build 
    esac
done

# get MTD folder place; same as Install.sh script file path (in the MTD folder)
MTDIR=$(dirname $(readlink -f $0))
cd $MTDIR # MTD folder place

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
kraken2-build --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --download-library archaea --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --download-library bacteria --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --download-library protozoa --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --download-library fungi --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --download-library plasmid --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --download-library UniVec_Core --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --add-to-library viruses4kraken.fa --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

# Bracken database update/re-building
if [[ $kmer == "" ]]; then
    bracken-build -d $MTDIR/kraken2DB_micro -t $threads -l $read_len
else
    bracken-build -d $MTDIR/kraken2DB_micro -t $threads -l $read_len -k $kmer
fi

echo "update done"