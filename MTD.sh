#!/bin/bash

while getopts i:o:h:t:c: option
do
    case "${option}" in
        i) inputdr=${OPTARG};;
        o) outputdr=${OPTARG};;
        h) hostid=${OPTARG};;
        t) threads=${OPTARG};;
    esac
done

# inputdr=~/RNAseq_raw_data/samplesheet.csv # select input directory; must store paried .fq.gz (eg. DJ01_1.fq.gz and DJ01_2.fq.gz) of each sample in the same folder as the samplesheet.csv
# outputdr=~/MTD_Results/test1 # select outputdr directory
# hostid=9544 # Enter host species taxonomy ID; initally supporting 9544 (rhesus monkey), 9606 (human), and 10090 (mouse).
# threads=20 # CPU threads; suggest >=16, eg. 20

# get MTD.sh script file path (in the MTD folder)
MTDIR=$(dirname $(readlink -f $0))
#parentname="$(dirname "$MTDIR")"
echo "MTD directory is $MTDIR"

# get conda path
condapath=$(head -n 1 $MTDIR/condaPath)
# activate MTD conda environment
source $condapath/etc/profile.d/conda.sh
conda activate MTD

inputdr=$(dirname $inputdr)
mkdir -p $outputdr
mkdir -p $outputdr/temp
cd $outputdr/temp

# Step 0: Host database auto selection
if [ $hostid == 9606 ]; then
    DB_host=$MTDIR/kraken2DB_human # for kraken2
    DB_hisat2=$MTDIR/hisat2_index_human/genome_tran #for hisat2
    gtf=$MTDIR/ref_human/Homo_sapiens.GRCh38.104.gtf.gz # for featureCounts

elif [ $hostid == 9544 ]; then
    DB_host=$MTDIR/kraken2DB_rhesus # for kraken2
    DB_hisat2=$MTDIR/hisat2_index_rhesus/genome_tran #for hisat2
    gtf=$MTDIR/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz # for featureCounts

elif [ $hostid == 10090 ]; then
    DB_host=$MTDIR/kraken2DB_mice
    DB_hisat2=$MTDIR/hisat2_index_mouse/genome_tran
    gtf=$MTDIR/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz # for featureCounts

elif [ -d "$MTDIR/kraken2DB_${hostid}" ]; then # test if customized host species exist
    DB_host=$MTDIR/kraken2DB_${hostid}
    DB_hisat2=$MTDIR/hisat2_index_${hostid}/genome_tran
    gtf=$MTDIR/ref_${hostid}/ref_${hostid}.gtf.gz

else
    echo "Host species is not supported. You can use bash Customized_host.sh for building."
    exit
fi
DB_micro=$MTDIR/kraken2DB_micro # customized kraken database for microbiome

# To extract sample names from input fastq files (support .fq.gz, .fastq.gz, .fq, or .fastq)
files=$(find $inputdr -name "*.fq.gz" -or -name "*.fastq.gz" -or -name "*.fq" -or -name "*.fastq" -type f)
files1=$(find $inputdr -name "*_1.fq.gz" -or -name "*_1.fastq.gz" -or -name "*_1.fq" -or -name "*_1.fastq" -type f)
#b=$(basename -a $inputdr) # store basenames of input directories into variable b to make a list of input sample names (eg. DJ01 EM77...)
for i in $files1; do
    fn=$(basename $i) #Extract file name, eg. DJ01_1.fq.gz
    sn=$(echo $fn | awk -F '_1' '{print $(NF-1)}') #Extract sample name, eg. DJ01
    lsn=$lsn" "$sn #Make a list of sample names; store basenames of input directories into variable lsn to make a list of input sample names (eg. DJ01 EM77...)
done

# Raw reads trimming
for i in $lsn; do # store input sample name in i; eg. DJ01
    # To get the corresponding fastq file as input (support .fq.gz, .fastq.gz, .fq, or .fastq)
    fq1=$(find $inputdr -type f \( -name "${i}_1.fq.gz" -or -name "${i}_1.fastq.gz" -or -name "${i}_1.fq" -or -name "${i}_1.fastq" \))
    fq2=$(find $inputdr -type f \( -name "${i}_2.fq.gz" -or -name "${i}_2.fastq.gz" -or -name "${i}_2.fq" -or -name "${i}_2.fastq" \))
	    #fastp with polyA/T trimming
        fastp --trim_poly_x \
            --length_required 40 \
            --thread 16 \
            -i $fq1 -I $fq2 \
            -o Trimmed_${i}_1.fq.gz -O Trimmed_${i}_2.fq.gz 
done

# Reads classification by kraken2; 1st step for host
for i in $lsn; do # store input sample name in i; eg. DJ01
    kraken2 --db $DB_host --use-names \
        --report Report_host_$i.txt \
        --threads $threads \
        --gzip-compressed \
        --paired \
        --classified-out ${i}_host#.fq \
        --unclassified-out ${i}_non-host_raw#.fq \
        Trimmed_${i}_1.fq.gz Trimmed_${i}_2.fq.gz \
        > Report_host_$i.kraken
done

# Reads classification by kraken2; 2nd step for non-host reads
for i in $lsn; do # store input sample name in i; eg. DJ01
    kraken2 --db $DB_micro --use-names \
        --report Report_non-host.raw_$i.txt \
        --threads $threads \
        --paired \
        --classified-out ${i}_raw_cseqs#.fq \
        --unclassified-out ${i}_raw_ucseqs#.fq \
        ${i}_non-host_raw_1.fq ${i}_non-host_raw_2.fq \
        > Report_non-host_raw_$i.kraken
done

# Decontamination step
conta_file=$MTDIR/conta_ls.txt
if test -f "$conta_file"; then
    tls=$(awk -F '\t' '{print $2}' $conta_file)
    conta_ls="${tls//$'\r\n'/ }"
    for i in $lsn; do
        python $MTDIR/Tools/KrakenTools/extract_kraken_reads.py \
            -k Report_non-host_raw_${i}.kraken \
            -s1 ${i}_non-host_raw_1.fq -s2 ${i}_non-host_raw_2.fq \
            -o ${i}_non-host_1.fq -o2 ${i}_non-host_2.fq \
            -r Report_non-host.raw_${i}.txt \
            --taxid $conta_ls --exclude --include-children
    done

    # Reads classification by kraken2; 3rd step for decontaminated non-host reads to get reports
    for i in $lsn; do
        kraken2 --db $DB_micro --use-names \
            --report Report_non-host_$i.txt \
            --threads $threads \
            --paired \
            --classified-out ${i}_cseqs#.fq \
            --unclassified-out ${i}_ucseqs#.fq \
            ${i}_non-host_1.fq ${i}_non-host_2.fq \
            > Report_non-host_$i.kraken
    done
fi

# Bracken analysis
for i in $lsn; do # store input sample name in i; eg. DJ01
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.phylum.bracken -r 75 -l P -t $threads
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.genus.bracken -r 75 -l G -t $threads
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.species.bracken -r 75 -l S -t $threads
done

#combined .bracken files (table like) into a single outputdr for Deseq2
python $MTDIR/Tools/combine_bracken_outputs.py --files *.phylum.bracken -o $outputdr/bracken_phylum_all
python $MTDIR/Tools/combine_bracken_outputs.py --files *.genus.bracken -o $outputdr/bracken_genus_all
python $MTDIR/Tools/combine_bracken_outputs.py --files *.species.bracken -o $outputdr/bracken_species_all

# move _bracken report files (tree like) to a separate folder
mkdir -p Report_non-host_bracken_species_normalized
mv *_bracken_species.txt Report_non-host_bracken_species_normalized
cd Report_non-host_bracken_species_normalized

#trim the name of _bracken report files (tree like) to the sample name (eg. DJ01)
for i in $lsn; do
    mv *${i}* $i
done

# Adjust bracken file (tree like) by normalizated reads counts; for additional visualization (.biom, .mpa, .krona)
Rscript $MTDIR/Normalization_afbr.R $outputdr/bracken_species_all $inputdr/samplesheet.csv $outputdr/temp/Report_non-host_bracken_species_normalized

#Converted _bracken report files (tree like) into .biom file for diversity analysis in phyloseq (R)
kraken-biom * -o $outputdr/bracken_species_all.biom --fmt json
#kraken-biom *_bracken_phylum -o bracken_phylum_all.biom --fmt json
#kraken-biom *_bracken_genus -o bracken_genus_all.biom --fmt json

# remove "sp. " in the .biom file; correct improper format before run export2graphlan.py
sed -i 's/sp. //g' $outputdr/bracken_species_all.biom

# go to temp folder
cd ../

mkdir -p ../graphlan
cd ../graphlan

#source $condapath/etc/profile.d/conda.sh
conda deactivate
conda activate py2

python $MTDIR/Tools/export2graphlan/export2graphlan.py \
    -i ../bracken_species_all.biom \
    -a annot.txt -t tree.txt \
    --discard_otus --most_abundant 50 \
    --annotations 2,3,4,5,6 \
    --external_annotations 7 --internal_levels --max_clade_size 300

conda deactivate
conda activate MTD

cd ../temp

# DEG & Annotation & Plots & Diversity & Preprocess
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/bracken_species_all $inputdr/samplesheet.csv $hostid

cd $outputdr/temp
mkdir -p bracken_raw_results # save the raw output from bracken (table like)
mv ../bracken_*_all bracken_raw_results

cd ../graphlan

python $MTDIR/Tools/graphlan/graphlan_annotate.py --annot annot.txt tree.txt outtree.txt # attach annotation to the tree
python $MTDIR/Tools/graphlan/graphlan.py --dpi 300 --size 7.0 outtree.txt outimg.png # generate the graphlan image

cd ../temp

## Visualization preprocess
# For krona
mkdir -p ../krona
for i in $lsn; do # store input sample name in i; eg. DJ01
    python $MTDIR/Tools/KrakenTools/kreport2krona.py \
        -r Report_non-host_bracken_species_normalized/${i} \
        -o ../krona/${i}-bracken.krona
done

# To make MPA style file
for i in $lsn; do # store input sample name in i; eg. DJ01
    python $MTDIR/Tools/KrakenTools/kreport2mpa.py \
        --display-header \
        -r Report_non-host_bracken_species_normalized/${i} \
        -o ${i}-bracken.mpa.txt
done
# Combine MPA files
python $MTDIR/Tools/KrakenTools/combine_mpa.py \
    -i *.mpa.txt \
    -o ../Combined.mpa

# HUMAnN3
mkdir -p HUMAnN_output
# HUMAnN does not run on paired-end reads by default; Concatenate fq files first
for n1 in *\_non-host_1.fq; do
    n2=${n1/_non-host_1/_non-host_2}
    cat $n1 $n2 > HUMAnN_output/$n1
done

cd HUMAnN_output

for file in *; do #trim the file name
    mv $file ${file/_non-host_1/}
done

# Run HUMAnN3
for i in *.fq; do
    humann --input $i \
        --output hmn_output \
        --threads $threads
done

#Join all gene family and pathway abudance files
humann_join_tables -i hmn_output/ -o humann_pathabundance.tsv --file_name pathabundance
humann_join_tables -i hmn_output/ -o humann_genefamilies.tsv --file_name genefamilies

# #Normalizing RPKs to CPM
# humann_renorm_table --input humann_pathabundance.tsv --output humann_pathabundance_cpm.tsv --units cpm --update-snames
# humann_renorm_table --input humann_genefamilies.tsv --output humann_genefamilies_cpm.tsv --units cpm --update-snames

#Normalizing RPKs to "relab" (relative abundance)
humann_renorm_table --input humann_pathabundance.tsv --output humann_pathabundance_relab.tsv --units relab --update-snames
humann_renorm_table --input humann_genefamilies.tsv --output humann_genefamilies_relab.tsv --units relab --update-snames

#Generate stratified tables; This utility will split a table into two files (one stratified and one unstratified).
humann_split_stratified_table --input humann_pathabundance_relab.tsv --output ./
humann_split_stratified_table --input humann_genefamilies_relab.tsv --output ./
    #Stratify unnormalized table (for Deseq2)
    humann_split_stratified_table --input humann_pathabundance.tsv --output ./
    humann_split_stratified_table --input humann_genefamilies.tsv --output ./

#Regroup gene familites table into KEGG orthologs and GO terms
humann_regroup_table --input humann_genefamilies_relab_stratified.tsv --groups uniref90_ko \
    --output humann_genefamilies_relAbundance_kegg.tsv
humann_regroup_table --input humann_genefamilies_relab_stratified.tsv --groups uniref90_go \
    --output humann_genefamilies_relAbundance_go.tsv
    #Regroup unnormalized table (for Deseq2)
    humann_regroup_table --input humann_genefamilies_stratified.tsv --groups uniref90_ko \
        --output humann_genefamilies_Abundance_kegg.tsv
    humann_regroup_table --input humann_genefamilies_stratified.tsv --groups uniref90_go \
        --output humann_genefamilies_Abundance_go.tsv

# Translate KEGG and GO ID to human readable terms
Rscript $MTDIR/humann_ID_translation.R \
    $outputdr/temp/HUMAnN_output/humann_genefamilies_relAbundance_kegg.tsv \
    $outputdr/temp/HUMAnN_output/humann_genefamilies_relAbundance_go.tsv \
    $MTDIR
    # Tranlate unnormalized table (for Deseq2)
    Rscript $MTDIR/humann_ID_translation.R \
        $outputdr/temp/HUMAnN_output/humann_genefamilies_Abundance_kegg.tsv \
        $outputdr/temp/HUMAnN_output/humann_genefamilies_Abundance_go.tsv \
        $MTDIR

#Cleaning up file structure
mkdir $outputdr/hmn_pathway_abundance_files
mkdir $outputdr/hmn_genefamily_abundance_files
mv *pathabundance* $outputdr/hmn_pathway_abundance_files/
mv *genefamilies* $outputdr/hmn_genefamily_abundance_files/

# #Translate KEGG and GO ID to human readable terms
# Rscript $MTDIR/humann_ID_translation.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_kegg.tsv \
#     $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_go.tsv

# DEG & Annotation & Plots & Diversity & Preprocess
cd $outputdr/hmn_genefamily_abundance_files
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_kegg_translated.tsv $inputdr/samplesheet.csv
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_go_translated.tsv $inputdr/samplesheet.csv

#humann_barplot
# humann_barplot --input $outputdr/hmn_pathway_abundance_files/humann_pathabundance_cpm_stratified.tsv \
#     --focal-metadatum Group --last-metadatum Group \
#     --focal-feature PWY-3781 \
#     --output $outputdr/hmn_pathway_abundance_files/humann_pathabundance_barplot.png
# humann_barplot --input $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_cpm_stratified.tsv \
#     --output $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_barplot.png

## continue to process the host reads
cd $outputdr/temp
# hisat2 alignment
for i in $lsn; do # store input sample name in i; eg. DJ01
    hisat2 -p $threads -q \
        -x $DB_hisat2 \
        --summary-file ${i}_hisat2_summary.txt \
        -1 ${i}_host_1.fq \
        -2 ${i}_host_2.fq \
        -S $i.sam
done

# featureCounts
featureCounts -T $threads \
   -p \
   -a $gtf \
   -o $outputdr/host_counts.txt \
   *.sam

for i in $lsn; do
    samtools view -bS $i.sam > $i.bam -@ $threads
    samtools sort $i.bam -o $i.sorted.bam -@ $threads
    samtools index $i.sorted.bam -@ $threads
done

mkdir -p BAM
mv *.sorted.bam *.sorted.bam.bai BAM/

cd $outputdr
# trim the featureCounts output(host_counts.txt) for downstream analysis
# delete the first line/row of a file then trim the sample name
sed '1d; 2 s/\.sam//g' host_counts.txt > tmpfile; mv tmpfile host_counts.txt

# DEG & Annotation & Plots & preprocess for host
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/host_counts.txt $inputdr/samplesheet.csv $hostid

# ssGSEA
Rscript $MTDIR/gct_making.R $outputdr/Host_DEG/host_counts_DEG.csv $inputdr/samplesheet.csv

Rscript $MTDIR/Tools/ssGSEA2.0/ssgsea-cli.R \
    -i $outputdr/ssGSEA/host.gct \
    -o $outputdr/ssGSEA/ssgsea-results \
    -d $MTDIR/Tools/ssGSEA2.0/db/msigdb/c2.all.v7.0.symbols.gmt \
    -y $MTDIR/Tools/ssGSEA2.0/config.yaml \
    -u $threads

Rscript $MTDIR/for_halla.R $outputdr/ssGSEA/ssgsea-results-scores.gct $inputdr/samplesheet.csv

echo "MTD DEG analyses are done. Starting microbiome x host association analyses..."

# halla: association analysis
#mkdir -p $outputdr/Associations
conda deactivate
conda activate halla0818
# for microbiome x host_genes
mkdir -p $outputdr/halla/host_gene # need to create a new directory for output to avoid "exists; deleting..." issue by halla
halla -x $outputdr/halla/Microbiomes.txt \
    -y $outputdr/halla/Host_gene.txt \
    -o $outputdr/halla/host_gene \
    --x_dataset_label Microbiomes \
    --y_dataset_label Host_gene \
    --diagnostic_plot -m spearman

# show all clusters
hallagram \
    -i $outputdr/halla/host_gene \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label Microbiomes \
    --y_dataset_label Host_gene \
    --output $outputdr/halla/host_gene/hallagram_all.png \
    --block_num -1

# for microbiome x host_pathways(ssGSEA)
mkdir -p $outputdr/halla/pathway
halla -x $outputdr/halla/Microbiomes.txt \
    -y $outputdr/halla/Host_score.txt \
    -o $outputdr/halla/pathway \
    --x_dataset_label Microbiomes \
    --y_dataset_label Host_pathway \
    --diagnostic_plot -m spearman

# show all clusters
hallagram \
    -i $outputdr/halla/pathway \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label Microbiomes \
    --y_dataset_label Host_pathway \
    --output $outputdr/halla/pathway_hallagram_all.png \
    --block_num -1

echo "MTD running is finished"