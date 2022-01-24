# MTD: Meta-Transcriptome Detector
MTD is a software that has two sub-pipelines to jointly analyze of host transcriptome with its microbiome by using bulk RNA-seq and single-cell RNA-seq data, respectively. It supports comprehensive microbiome species and vectors, including viruses, bacteria, protozoa, fungi, plasmids, and vectors. MTD is executed in Bash in GNU/Linux system. Users can easily install and run MTD using only one command and without requiring root privileges. The outputs (graphs, tables, count matrixes, etc.) are automatically generated and stored in the designated directory/folder defined by the user.
# Key Points
* MTD enables **simultaneous analyses** of the microbiome and the host cell transcriptome in **bulk** and **single-cell RNA-seq data**.
* The **correlation** between the microbiome and the host transcriptome can be automatically analyzed.
* MTD has an extensive microbiome detection capacity, including **viruses, bacteria, protozoa, fungi, plasmids, and vectors**.
* Installation and use MTD is as easy as **one command line** without the requirement of administrator/root privilege.
* **Decontamination function** is enabled to eliminate the common contaminant microbes in the laboratory environment.
# Requirements
* 160 Gb RAM
* 560 Gb storage space
* Conda was installed
* GNU/Linux system with Bash
# Installation
1. Download directly or git clone MTD (git clone https://github.com/FEI38750/MTD.git) to the place you want to install. The full version of software (~530Gb) will be installed in this MTD folder.
2. In termial, type\
**bash [path/to/MTD]/Install.sh -t [threads] -p [path/to/conda]**\
For example:
        <pre><code>bash ~/MTD/Install.sh -t 20 -p ~/miniconda3
</code></pre>
## Notes
* Installation may take 1-2 days (10-20 threads).
* If conda hasn't been installed in your system, please use the code below to install the conda:
<addr><pre><code>URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
curl $URL > miniconda-installer.sh
bash miniconda-installer.sh -b</code></pre>
Conda then will be installed in your home directory, such as path: ~/miniconda3
* MTD conda environments occupies ~24Gb in your conda folder.
* Tips: file management software such as FileZilla (https://filezilla-project.org/download.php?show_all=1) can help you to manage your files on HPC/server.
# Run MTD
## Bulk RNA-seq
### Preparation of the raw data (local data)
1. Store all the paired-end fastq files (accepted: fastq, fastq.gz, fq, fq.gz) to be analyzed in a folder, subfolders for each sample are accepted.
  The paired fastq files must be named starting with the sample name followed by "_1" and "_2". For example, sample1_1.fq.gz and sample1_2.fq.gz are paired-end fastq files for sample1.
&nbsp;&nbsp;<img src="https://github.com/FEI38750/MTD/blob/main/Img/input_folder1.jpg">
2. Prepare the samplesheet.csv. You can copy and modify the one in MTD folder.
  ![image1](https://github.com/FEI38750/MTD/blob/main/Img/Tutorial1.jpg)
3. Put samplesheet.csv in the same folder as the fastq files.\
&nbsp;<img src="https://github.com/FEI38750/MTD/blob/main/Img/input_folder.jpg" width=50%>
### Preparation of the raw data (online SRA data)
Alternatively, MTD can directly use data from SRA NCBI as input samples. Users just need to enter the corresponding SRR accessions into the sample_name column of samplesheet.csv \
<img src="https://github.com/FEI38750/MTD/blob/main/Img/SRR_bulk.png" width=60%> \
SRR inputs samples from NCBI can be downloaded and prepared automatically.
### Run
In termial, type\
  **bash [path/to/MTD]/MTD.sh -i [path/to/samplesheet.csv] -o [path/to/output_folder] -h [host species taxonomy ID] -t [threads]**\
Host species taxonomy ID: human:9606, mouse:10090, rhesus monkey:9544\
For example:
        <pre><code>bash ~/MTD/MTD.sh -i ~/raw_data/samplesheet.csv -o ~/MTD_output -h 9544 -t 20</code></pre>
### Notes
* Test run: For user who does not have the bulk RNA-seq raw data on hand could have a test run by command:
  <pre><code>bash [path/to/MTD]/MTD.sh -i [path/to/MTD]/test/Bulk_RNAseq/samplesheet.csv -o [path/to/MTD]/test/Bulk_RNAseq/output -h 9606 -t [threads]</code></pre>
* Users who prefer using Magic-BLAST instead of HISAT2 for host reads mapping can add <code>-b blast</code> flag. For example:
  <pre><code>bash ~/MTD/MTD.sh -i ~/raw_data/samplesheet.csv -o ~/MTD_output -h 9544 -t 20 -b blast</code></pre>
  
## Single-cell RNA-seq
1. Put the count matrix of host genes in a folder named with the sample name. In this folder, 10x should be a matrix.mtx, a genes.tsv, and a bar
        s.tsv; or a single .h5 file. Dropseq should be a .dge.txt file.
&nbsp;&nbsp;<img src="https://github.com/FEI38750/MTD/blob/main/Img/input_folder_SC1.jpg" width=70%>
2. Type the path of the host matrix folder and the corresponding fastq files into the columns of the samplesheet_SC.csv accordingly. For example:
&nbsp;&nbsp;<img src="https://github.com/FEI38750/MTD/blob/main/Img/samplesheet_SC1.jpg">
  Then MTD will read the corresponding file paths from this samplesheet_SC.csv for single-cell analysis.
3. In termial, type\
  **bash [path/to/MTD]/MTD_singleCell.sh -i [path/to/samplesheet_SC.csv] -o [path/to/Output_folder] -h [Host species taxonomy ID] -t [Threads] -p [Platform] -d [prime Direction] -c [path/to/Cell_barcode_file.whitelist.txt]**\
  Single cell RNAseq platform(-p): enter 1 for 10x or 2 for Dropseq platform\
  prime_direction(-d): specifying barcode locations: enter 3 or 5 for barcodes are at the 3’ end or 5' end of the read\
  For example:
        <pre><code>bash ~/MTD/MTD_singleCell.sh -i ~/scRNAseq_rawData/samplesheet_SC.csv -o ~/output -h 10090 -t 20 -p 1 -d 3</code></pre>
### Notes
* 10x and Dropseq use paired end sequence. The first fastq file contains barcodes (e.g., 26bp length in SRR4210_R1.fastq). The second fastq file contains transcript's sequences (e.g., 98bp length in SRR4210_R2.fastq).
* Default QC is *subset= nFeature_RNA>200 & nFeature_RNA < 2\*median(number_of_Feature_RNA) & percent.mt < 10*\
  In addition, user can customize QC by adding -l [Minimum nFeature_RNA] -r [Maximum nFeature_RNA] -m [percent.mt]
  
# Outputs
  ## Bulk RNA-seq
  ![image1](https://github.com/FEI38750/MTD/blob/main/Img/MTD_bulk.png)
  The results are generated automatically and saved in the output folder defined by the user.\
  The output included:\
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/Output_folder.jpg" width=70% height=70%>
* For **host**: [path/to/output_folder]/Host_DEG/\
The count matrix (host_counts_DEG.csv) contains the Ensembl gene ID, gene symbol, chromosome name, gene position, functional descriptions, DEG results for each group comparison, raw read counts, normalized reads count, normalized and transformed reads counts. This comprehensive count matrix facilitates the user to perform downstream analyses such as pathway enrichment and customized data visualization.\
The data visualization includes the heatmap (with/without gene name), Venn Diagram, PCA, barplot, and volcano plots.\
The individual group comparison results are saved in the corresponding subfolder (e.g., group1_vs_group2).\
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/bulk_output_folder.jpg" width=100% height=100%>
* For **microbiome**: [path/to/output_folder]/Nonhost_DEG/\
  The count matrix (bracken_normalized_species_all_DEG.csv) contains the name and taxonomy ID of microbiome species, DEG results for each group comparison, raw read counts, normalized reads count, normalized and transformed reads counts.\
Diversity analysis, unclassified reads comparison, abundance&DEG heatmaps, phylogenetic trees.\
Venn Diagram, heatmap, PCoA, barplot, and volcano plots for the results of species abundance and group comparisons.\
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/Nonhost_output_folder.jpg" width=70% height=70%>
* **Microbiome metabolic molecules**: hmn_genefamily_abundance_files contain microbiome metabolic molecules and group comparison results. Results are translated to kegg and go terms to facilitate reading and demonstrated via Venn Diagram, heatmap, PCA, barplot, and volcano plots and count matrix.\
hmn_pathway_abundance_files contain pathway results of those molecules (e.g. humann_pathabundance_relab_stratified.tsv contains normalized relative abundance of pathways).
* **Association analysis**: halla folder contains the results of association between:\
host gene and microbiome species\
host pathways and microbiome species
 ## Single-cell RNA-seq
  **Count matrix** for the single-cell microbiome is automatically generated and saved in the output folder.\
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/MTD_SC_countMatrix_demo.jpg" width=85% height=85%>\
  The results of correlation test between microbiome and host genes are generated automatically and saved in the output folder defined by the user.
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/SC_corr.jpg">
  
## Notes
  - For reference, the MTD running time would be:
    - Bulk RNA-seq: for 10 samples, 20 fastq files (total 47 Gb) by using 20 threads CPU is ~8 hours (except correlation analysis). In addition, the correlation analysis may need a further ~26-30hours (results in halla folder). So the total running time would be ~34-38 hours.
    - Single-cell RNA-seq: for fastq files contain ~2000 cells, by using 20 threads CPU is ~2 hours.
  - In addition to initially supported human, mouse, and rhesus monkey, users can easily add any other host species by one command line:\
    **bash Customized_host.sh -t [threads] -d [host_genome_Ensembl_address] -c [host_taxid] -g [host_gtf_Ensembl_address]**\
      For example:\
      <code>bash ~/MTD/Customized_host.sh -t 20 -d http://ftp.ensembl.org/pub/release-104/fasta/callithrix_jacchus/dna/Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz -c 9483 -g http://ftp.ensembl.org/pub/release-104/gtf/callithrix_jacchus/Callithrix_jacchus.ASM275486v1.104.gtf.gz </code>
  - MTD supported ~200 host species. Please refer to HostSpecies.csv for their Taxonomy IDs and names.
  - Users can update microbiome databases easily by one command line: **bash Update.sh -t [threads]**\
      For example:\
    <code>bash ~/MTD/Update.sh -t 20</code>
  - Users can modify the contaminant list (conta_ls.txt) in the MTD folder by adding the taxonomy ID of the microbe in the second column of the list and its name in the first column (optional).
## Advanced options
  - Users can provide additional metadata information of the samples (metadata.csv) for more complex experimental designs and analysis for bulk RNAseq data. In addition to sample_name and group, users can add the additional columns that represent covariants for the analysis. So the covariants, such as batch, group, age, sex, etc., can be adjusted during the read count abundance analysis through DESeq2, ANCOM-BC, and MaAsLin2 in MTD. Further information about a proper experimental design for analysis, please refer to:
https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html \
Prepare the metadata.csv. You can copy and modify the one in MTD folder.\
        <img src="https://github.com/FEI38750/MTD/blob/main/Img/metadata.jpg" width=60% height=60%> \
Then just to add the flag <code>-m [path/to/metadata.csv]</code> in the MTD bulk RNAseq command line. Such as:\
        <code>bash [path/to/MTD]/MTD.sh -i [path/to/samplesheet.csv] -o [path/to/output_folder] -h [host species taxonomy ID] -t [threads] -m [path/to/metadata.csv]</code>
  - For users who want to tune parameters, MTD has the optional settings by additional flags for the potential important steps. Please refer to the MTD/Tutorial/Advanced_options.xlsx for explanations.
  - For users who have advanced knowledge and wants to have further complicated settings, could add or change corresponding parameters inside the source code MTD/MTD.sh, for example, by searching *# fastp* to locate the code block of fastp settings then add additional parameters according to the options on https://github.com/OpenGene/fastp#all-options. Searching *# HISAT2* or *# Magic-BLAST* to the place for parameters adjustion for the host reads mapping.
  - For users who run the MTD job on HPC interactively, it is optimal running through the Linux GNU Screen tool, which can prevent the interruption due to the user end (e.g. internet disconnect). Screen has already installed in the MTD Conda environment, user first type: <code>bash conda activate MTD</code> to activate the environment.\
    To start a screen session, simply type screen in your console: <code>screen</code> \
    To resume your screen session use the following command: <code>screen -r</code> \
    To close your screen session use: <code>screen -X -S [session # you want to kill] quit</code>
  
# Citation
Fei Wu, Yao-Zhong Liu, Binhua Ling. (2021). MTD: a unique pipeline for host and meta-transcriptome joint and integrative analyses of RNA-seq data. bioRxiv 2021.11.16.468881; doi: https://doi.org/10.1101/2021.11.16.468881
  
# Licence
This software is freely available for academic users. Usage for commercial purposes is not allowed. Please refer to the LICENCE page.
