#!/bin/bash

# default settings
PD=""
scn=""
length=40 # read length trimming by fastp

while getopts i:o:h:t:p:d:l:r:m:f: option
do
    case "${option}" in
        i) samplesheet=${OPTARG};;
        o) outputdr=${OPTARG};;
        h) hostid=${OPTARG};;
        t) threads=${OPTARG};;
        p) platform=${OPTARG};;
        d) primeD=${OPTARG};;
        l) low_RNA=${OPTARG};;
        r) high_RNA=${OPTARG};;
        m) mt_percent=${OPTARG};;
        f) length=${OPTARG};;
    esac
done
echo "input samlesheet_SC.csv: $samplesheet"
echo "output directory: $outputdr"
echo "host taxonomy id: $hostid"
echo "CPU threads: $threads"
echo "single cell sequencing platform: $platform"
echo "barcode locations(on 5' or 3' end of reads): $primeD"

# inputdr=~/scRNAseq/GSE121611/samplesheet_SC.csv # select input path for samplesheet.csv
# outputdr=~/scRNAseq/GSE121611/output # select output directory
# threads=20 # CPU threads; suggest >=16, eg. 20
# hostid=10090 # Enter host species taxonomy ID; MTD supports 9606/9544/10090 initially
# platform=1 # enter 1 for 10x or 2 for dropseq platform
# prime_direction=3 # specifying barcode locations: enter 3 or 5 for barcodes are at the 3’ end or 5' end of the read

mkdir -p $outputdr
mkdir -p $outputdr/temp
cd $outputdr/temp

if [[ $primeD == 3 ]]; then
    PD="--3prime "
    echo "prime_direction is $PD"
else
    echo "prime_direction is 5prime(default)"
fi

# get MTD.sh script file path (in the MTD folder)
MTDIR=$(dirname $(readlink -f $0))
#parentname="$(dirname "$MTDIR")"
echo "MTD directory is $MTDIR"

# get conda path
condapath=$(head -n 1 $MTDIR/condaPath)
# activate MTD conda environment
source $condapath/etc/profile.d/conda.sh

conda activate R412
echo 'Single cell fastq pre-processing: arrangement and whilelist making...'
Rscript $MTDIR/SingleCell_Prep.R $samplesheet $outputdr $platform
conda deactivate

conda activate MTD
echo 'MTD running  progress:'
echo '>>                  [10%]'

inputdr=$outputdr/fastq

# Step 0: Host database auto selection
if [[ $hostid == 9606 ]]; then
    DB_host=$MTDIR/kraken2DB_human # for kraken2
elif [[ $hostid == 9544 ]]; then
    DB_host=$MTDIR/kraken2DB_rhesus # for kraken2
elif [[ $hostid == 10090 ]]; then
    DB_host=$MTDIR/kraken2DB_mice # for kraken2
elif [[ -d "$MTDIR/kraken2DB_${hostid}" ]]; then # test if customized host species exist
    DB_host=$MTDIR/kraken2DB_${hostid}
else
    echo "Host species is not supported. Please use bash Customized_host.sh to build one."
    exit
fi

DB_micro=$MTDIR/kraken2DB_micro # customized kraken database for microbiome

# To set bar code pattern
if [[ $platform == 1 ]]; then
    BCpattern=CCCCCCCCCCCCCCCCNNNNNNNNNN
elif [[ $platform == 2 ]]; then
    BCpattern=CCCCCCCCCCCCNNNNNNNN
else
    echo "Single cell sequencing platform is not supported"
    exit
fi

# To extract sample names from input fastq files (support .fq.gz, .fastq.gz, .fq, or .fastq)
files1=$(find $inputdr -name "*_1.fq.gz" -or -name "*_1.fastq.gz" -or -name "*_1.fq" -or -name "*_1.fastq" -type f)
files2=$(find $inputdr -name "*_2.fq.gz" -or -name "*_2.fastq.gz" -or -name "*_2.fq" -or -name "*_2.fastq" -type f)
# To make a array(list) of sample names
lsn=()
for i in $files1; do
    fn=$(basename $i)
    sn=$(echo $fn | awk -F '_1' '{print $(NF-1)}')
    lsn+=( $sn )
done

echo 'Step 1: Extract barcdoes and UMIs and add to read names...'
for i in "${lsn[@]}"
do
    fq1=$(find $inputdr -type f \( -name "${i}_1.fq.gz" -or -name "${i}_1.fastq.gz" -or -name "${i}_1.fq" -or -name "${i}_1.fastq" \))
    fq2=$(find $inputdr -type f \( -name "${i}_2.fq.gz" -or -name "${i}_2.fastq.gz" -or -name "${i}_2.fq" -or -name "${i}_2.fastq" \))
    # Get the cell barcodes (whitelist)
    CB=$(find ${outputdr}/Cell_Barcode -type f -name "${i}_Barcodes.tsv" )
    # Extract barcdoes and UMIs and add to read names
    umi_tools extract --extract-method=string \
                    --bc-pattern=$BCpattern ${PD}\
                    --stdin $fq1 \
                    --stdout ${i}_1.extracted.fastq.gz \
                    --read2-in $fq2 \
                    --read2-out=${i}_2.extracted.fastq.gz \
                    --whitelist=$CB
done

echo 'MTD running  progress:'
echo '>>>>                [20%]'

echo 'Step 2: Raw reads trimming and quality control...'
for i in "${lsn[@]}" # store input sample name in i; eg. SRR7819592
do
    #fastp with polyA/T trimming
        fastp --trim_poly_x \
            --length_required $length \
            --thread $threads \
            -i ${i}_2.extracted.fastq.gz \
            -o ${i}_2.extracted.trimmed.fastq.gz
done

echo 'MTD running  progress:'
echo '>>>>>>              [30%]'

echo 'Step 3: Classify reads by kraken2...'
echo 'Reads classification by kraken2; 1st step for host classification...'
for i in "${lsn[@]}"; do
    kraken2 --db $DB_host --use-names \
        --report Report_host_$i.txt \
        --threads $threads \
        --gzip-compressed \
        --classified-out ${i}_host.fq \
        --unclassified-out ${i}_non-host_raw.fq \
        ${i}_2.extracted.trimmed.fastq.gz \
        > Report_host_$i.kraken
done

echo 'MTD running  progress:'
echo '>>>>>>>>            [40%]'

echo 'Reads classification by kraken2; 2nd step for non-host reads classification...'
for i in "${lsn[@]}"; do
    kraken2 --db $DB_micro --use-names \
        --report Report_non-host.raw_$i.txt \
        --threads $threads \
        --classified-out ${i}_raw_cseqs.fq \
        --unclassified-out ${i}_raw_ucseqs.fq \
        ${i}_non-host_raw.fq \
        > Report_non-host_raw_$i.kraken
done

echo 'MTD running  progress:'
echo '>>>>>>>>>>          [50%]'

echo 'Step 4: Decontamination...'
conta_file=$MTDIR/conta_ls.txt
if test -f "$conta_file"; then
    tls=$(awk -F '\t' '{print $2}' $conta_file)
    conta_ls="${tls//$'\r\n'/ }"
    for i in "${lsn[@]}"; do
        python $MTDIR/Tools/KrakenTools/extract_kraken_reads.py \
            -k Report_non-host_raw_${i}.kraken \
            -s1 ${i}_non-host_raw.fq \
            -o ${i}_non-host.fq \
            -r Report_non-host.raw_${i}.txt \
            --taxid $conta_ls --exclude --include-children
    done

    echo 'MTD running  progress:'
    echo '>>>>>>>>>>>>        [60%]'

    echo 'Reads classification by kraken2; 3rd step for decontaminated non-host reads to get reports...'
    for i in "${lsn[@]}"; do
        kraken2 --db $DB_micro --use-names \
            --report Report_non-host_$i.txt \
            --threads $threads \
            --classified-out ${i}_cseqs.fq \
            --unclassified-out ${i}_ucseqs.fq \
            ${i}_non-host.fq \
            > Report_non-host_$i.kraken
    done
fi

echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>      [70%]'

echo 'Step 5: Count UMIs per gene per cell for microbiome...'
for file in Report_non-host_*.kraken
do
    awk -F "_" 'gsub("\t","_") {print $2"_"$4"_"$3"\t"$5}' $file | sort -t $'\t' -k2 | \
    umi_tools count_tab --per-cell -S $file.c.tsv
done

echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>    [80%]'

conda deactivate
conda activate R412
echo 'Step 6: Make the count matrix of microbiome...'
Rscript $MTDIR/Singelcell4kraken2_umitools_batch.R $outputdr/temp

for i in "${lsn[@]}"; do
    mv Report_non-host_${i}.kraken.c.tsv_Count.txt $outputdr/${i}_count_matrix.txt
done

echo "Count matrix for microbiome is done"

echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>>>  [90%]'

echo 'Ex-step: Find out the microbiome and host gene with high correlation...'
Rscript $MTDIR/SC_corr.R $samplesheet $outputdr $platform $threads $low_RNA $high_RNA $mt_percent

echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>>>>>[100%]'
echo "Count matrix for microbiome is in the output folder with name end with _count_matrix.txt"
echo "MTD running is finished"