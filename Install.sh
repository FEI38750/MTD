#!/bin/bash

condapath=~/miniconda3
while getopts t:p: option # user can ingore -c and -p if conda has been already installed in the home directory
do
    case "${option}" in
        t) threads=${OPTARG};;
        p) condapath=${OPTARG};; # path to miniconda/anaconda (default is in home directory: ~/miniconda3)
    esac
done

# get MTD folder place; same as Install.sh script file path (in the MTD folder)
dir=$(dirname $(readlink -f $0))
cd $dir # MTD folder place

touch condaPath
echo "$condapath" > $dir/condaPath

source $condapath/etc/profile.d/conda.sh

conda env create -f Installation/MTD.yml
conda env create -f Installation/py2.yml
conda env create -f Installation/halla0818.yml

conda activate py2 # install dependencies of py2 in case pip does work in conda yml
pip install backports-functools-lru-cache==1.6.1 biom-format==2.0.1 cycler==0.10.0 h5py==2.10.0 hclust2==1.0.0 kiwisolver==1.1.0 matplotlib==2.2.5 numpy==1.16.6 pandas==0.24.2 pyparsing==2.4.7 pyqi==0.3.2 python-dateutil==2.8.1 pytz==2021.1 scipy==1.2.3 six==1.15.0 subprocess32==3.5.4
conda deactivate

conda activate halla0818 # install dependencies of halla
R -e 'install.packages(c("XICOR","mclust","BiocManager"), repos="http://cran.us.r-project.org")'
R -e 'BiocManager::install("preprocessCore", ask = FALSE)'
R -e 'install.packages("eva", INSTALL_opts = "--no-lock", repos="http://cran.us.r-project.org")'
conda deactivate

conda activate MTD
wget -c https://www.genome.jp/ftp/db/virushostdb/virushostdb.genomic.fna.gz
unpigz virushostdb.genomic.fna.gz
cat Installation/M33262_SIVMM239.fa virushostdb.genomic.fna > viruses4kraken.fa

# Kraken2 database building - Microbiome
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

# Kraken2 database building - Human
DBNAME=kraken2DB_human
kraken2-build --download-taxonomy --threads $threads --db $DBNAME
kraken2-build --download-library human --threads $threads --db $DBNAME
kraken2-build --build --threads $threads --db $DBNAME

# Kraken2 database building - Mouse
DBNAME=kraken2DB_mice
mkdir -p $DBNAME
cd $DBNAME
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
unpigz GCF_000001635.27_GRCm39_genomic.fna.gz
mv GCF_000001635.27_GRCm39_genomic.fna GCF_000001635.27_GRCm39_genomic.fa
cd ..
kraken2-build --download-taxonomy --threads $threads --db $DBNAME
kraken2-build --add-to-library $DBNAME/GCF_000001635.27_GRCm39_genomic.fa --threads $threads --db $DBNAME
kraken2-build --build --threads $threads --db $DBNAME

# Kraken2 database building - Rhesus macaque
DBNAME=kraken2DB_rhesus
mkdir -p $DBNAME
cd $DBNAME
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna.gz
unpigz GCF_003339765.1_Mmul_10_genomic.fna.gz
mv GCF_003339765.1_Mmul_10_genomic.fna GCF_003339765.1_Mmul_10_genomic.fa
cd ..
kraken2-build --download-taxonomy --threads $threads --db $DBNAME
kraken2-build --add-to-library $DBNAME/GCF_003339765.1_Mmul_10_genomic.fa --threads $threads --db $DBNAME
kraken2-build --build --threads $threads --db $DBNAME

# Bracken database building
bracken-build -d $dir/kraken2DB_micro -t $threads -l 75

# install HUMAnN3 databases
mkdir -p $dir/HUMAnN/ref_database/
cd $dir/HUMAnN/ref_database/
wget -c http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan.v296_201901.tar.gz
wget -c http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_annotated_v201901.tar.gz
wget -c http://huttenhower.sph.harvard.edu/humann2_data/full_mapping_v201901.tar.gz

mkdir -p $dir/HUMAnN/ref_database/chocophlan
tar xzvf full_chocophlan.v296_201901.tar.gz -C chocophlan
mkdir -p $dir/HUMAnN/ref_database/full_UniRef90
tar xzvf uniref90_annotated_v201901.tar.gz -C full_UniRef90
mkdir -p $dir/HUMAnN/ref_database/utility_mapping
tar xzvf full_mapping_v201901.tar.gz -C utility_mapping
cd $dir

humann_config --update database_folders nucleotide $dir/HUMAnN/ref_database/chocophlan
humann_config --update database_folders protein $dir/HUMAnN/ref_database/full_UniRef90
humann_config --update database_folders utility_mapping $dir/HUMAnN/ref_database/utility_mapping

# install host references
# download host GTF
    # download rhesus macaque GTF
    wget -c http://ftp.ensembl.org/pub/release-104/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.104.gtf.gz -P ref_rhesus
    # download human GTF
    wget -c http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz -P ref_human
    # download mouse GTF
    wget -c http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz -P ref_mouse

# Building indexes for hisat2
# rhesus macaques
mkdir -p hisat2_index_rhesus
cd hisat2_index_rhesus
cp ../ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz .
gzip -d Macaca_mulatta.Mmul_10.104.gtf.gz
mv Macaca_mulatta.Mmul_10.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
wget -c http://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz #use ensembl genome to compatible with featureCount
gzip -d Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
mv Macaca_mulatta.Mmul_10.dna.toplevel.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

# mouse
mkdir -p hisat2_index_mouse
cd hisat2_index_mouse
cp ../ref_mouse/Mus_musculus.GRCm39.104.gtf.gz .
gzip -d Mus_musculus.GRCm39.104.gtf.gz
mv Mus_musculus.GRCm39.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
wget -c http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz #use ensembl genome to compatible with featureCount
gzip -d Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
mv Mus_musculus.GRCm39.dna.primary_assembly.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

# human
mkdir -p hisat2_index_human
cd hisat2_index_human
cp ../ref_human/Homo_sapiens.GRCh38.104.gtf.gz .
gzip -d Homo_sapiens.GRCh38.104.gtf.gz
mv Homo_sapiens.GRCh38.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
wget -c http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz #use ensembl genome to compatible with featureCount
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

# # download preduild index for hisat2 # prebuild is from NCBI, may be not compatiable with featureCount
#     # H. sapiens
#     mkdir -p hisat2_index_human
#     cd hisat2_index_human
#     wget https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz
#     pigz -dc grch38_tran.tar.gz | tar xf -
#     cd ..

# install R packages
Rscript $dir/Installation/R_packages_installation.R

chmod +x MTD.sh

echo "MTD installation is finished"
