# MTD: Meta-Transcriptome Detector
MTD is a software that has two sub-pipelines to detect and quantify microbiomes by analyzing bulk RNA-seq data and single-cell RNA-seq data, respectively. It supports comprehensive microbiome species and vectors, including viruses, bacteria, protozoa, fungi, plasmids, and vectors. MTD is executed in Bash in GNU/Linux system. Users can easily install and run MTD using only one command and without requiring root privileges. The outputs (graphs, tables, count matrixes, etc.) are automatically generated and stored in the designated directory/folder defined by the user.
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
​       <pre><code>bash ~/MTD/Install.sh -t 20 -p ~/miniconda3
</code></pre>
## Notes
* Installation may take 1-2 days.
* If conda hasn't been installed in your system, please use the code below to install the conda:
<addr><pre><code>URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
curl $URL > miniconda-installer.sh
bash miniconda-installer.sh -b</code></pre>
Conda then will be installed in your home directory, such as path: ~/miniconda3
* MTD conda environments occupies ~24Gb in your conda folder.
* Tips: file management software such as FileZilla (https://filezilla-project.org/download.php?show_all=1) can help you to manage your files on HPC/server.
# Run MTD
## Bulk RNA-seq
1. Store all the paired-end fastq files (accepted: fastq, fastq.gz, fq, fq.gz) to be analyzed in a folder, subfolders for each sample are accepted.
  The paired fastq files must be named starting with the sample name followed by "_1" and "_2". For example, sample1_1.fq.gz and sample1_2.fq.gz are paired-end fastq files for sample1.
2. Prepare the samplesheet.csv. You can copy and modify the one in MTD folder.
  ![image1](https://github.com/FEI38750/MTD/blob/main/Img/Tutorial1.jpg)
3. Put samplesheet.csv in the same folder as the fastq files.
  ![image2](https://github.com/FEI38750/MTD/blob/main/Img/Tutorial2.jpg)
4. In termial, type\
  **bash [path/to/MTD]/MTD.sh -i [path/to/samplesheet.csv] -o [path/to/output_folder] -h [host species taxonomy ID] -t [threads]**\
Host species taxonomy ID: human:9606, mouse:10090, rhesus monkey:9544\
For example:
​       <pre><code>bash ~/MTD/MTD.sh -i ~/raw_data/samplesheet.csv -o ~/MTD_output -h 9544 -t 20</code></pre>
## Single-cell RNA-seq
1. Put the count matrix of host genes in a folder named with the sample name. In this folder, 10x should be a matrix.mtx, a genes.tsv, and a barcodes.tsv; or a single .h5 file. Dropseq should be a .dge.txt file.
2. Type this folder path into the column host_matrix_folder of the samplesheet_SC.csv. For example:
  ![imageSH](https://github.com/FEI38750/MTD/blob/main/Img/samplesheet_SC.jpg)
3. In termial, type\
  **bash [path/to/MTD]/MTD_singleCell.sh -i [path/to/Input_folder] -o [path/to/Output_folder] -h [Host species taxonomy ID] -t [Threads] -p [Platform] -d [prime Direction] -c [path/to/Cell_barcode_file.whitelist.txt]**\
  Single cell RNAseq platform(-p): enter 1 for 10x or 2 for Dropseq platform\
  prime_direction(-d): specifying barcode locations: enter 3 or 5 for barcodes are at the 3’ end or 5' end of the read\
  For example:
​       <pre><code>bash ~/MTD/MTD_singleCell.sh -i ~/scRNAseq_rawData/samplesheet_SC.csv -o ~/output -h 10090 -t 20 -p 1 -d 3</code></pre>
### Notes
* 10x and Dropseq use paired end sequence. The first fastq file contains barcodes (e.g., 26bp length in SRR4210_R1.fastq). The second fastq file contains transcript's sequences (e.g., 98bp length in SRR4210_R2.fastq).
* Default QC is *subset= nFeature_RNA>200 & nFeature_RNA < 2\*median(number_of_Feature_RNA) & percent.mt < 10*\
  In addition, user can customize QC by adding -l [Minimum nFeature_RNA] -r [Maximum nFeature_RNA] -m [percent.mt]
  
# Outputs
  ## Bulk RNA-seq
  ![image1](https://github.com/FEI38750/MTD/blob/main/Img/MTD_bulk.jpg)
  The results are generated automatically and saved in the output folder defined by the user.\
  The output included:\
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/Output_folder.jpg" width=70% height=70%>
* For **host**: [path/to/output_folder]/Host_DEG/\
The count matrix (host_counts_DEG.csv) contains the Ensembl gene ID, gene symbol, chromosome name, gene position, functional descriptions, DEG results for each group comparison, raw read counts, normalized reads count, normalized and transformed reads counts. This comprehensive count matrix facilitates the user to perform downstream analyses such as pathway enrichment and customized data visualization.\
The data visualization includes the heatmap (with/without gene name), Venn Diagram, PCA, barplot, and volcano plots.\
The individual group comparison results are saved in the corresponding subfolder (e.g., group1_vs_group2).\
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/Host_DEG.jpg" width=70% height=70%>
* For **microbiome**: [path/to/output_folder]/Nonhost_DEG/\
  The count matrix (bracken_normalized_species_all_DEG.csv) contains the name and taxonomy ID of microbiome species, DEG results for each group comparison, raw read counts, normalized reads count, normalized and transformed reads counts.\
Diversity analysis, unclassified reads comparison, abundance&DEG heatmaps, phylogenetic trees.\
Venn Diagram, heatmap, PCoA, barplot, and volcano plots for the results of species abundance and group comparisons.\
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/Nonhost_DEG.jpg" width=70% height=70%>
* **Microbiome metabolic molecules**: hmn_genefamily_abundance_files contain microbiome metabolic molecules and group comparison results. Results are translated to kegg and go terms to facilitate reading and demonstrated via Venn Diagram, heatmap, PCA, barplot, and volcano plots and count matrix.\
hmn_pathway_abundance_files contain pathway results of those molecules.
* **Association analysis**: halla folder contains the results of association between:\
host gene and microbiome species\
host pathways and microbiome species
 ## Single-cell RNA-seq
  **Count matrix** for the single-cell microbiome is automatically generated and saved in the output folder.\
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/SingleCell_results.jpg" width=80% height=80%>
  The results of correlation test between microbiome and host genes are generated automatically and saved in the output folder defined by the user.
  <img src="https://github.com/FEI38750/MTD/blob/main/Img/SC_corr.jpg">
  
## Notes
  - For reference, the MTD running time would be:
    - Bulk RNA-seq: for 10 samples, 20 fastq files (total 47 Gb) by using 20 threads CPU is ~8 hours (except correlation analysis). In addition, the correlation analysis may need a further ~26-30hours (results in halla folder). So the total running time would be ~34-38 hours.
    - Single-cell RNA-seq: for fastq files contain ~2000 cells, by using 20 threads CPU is ~2 hours.
  - Users can add any other host species easily by one command line:\
    **bash Customized_host.sh -t [threads] -d [host_genome_Ensembl_address] -c [host_taxid] -g [host_gtf_Ensembl_address]**\
      For example:\
      <code>bash ~/MTD/Customized_host.sh -t 20 -d http://ftp.ensembl.org/pub/release-104/fasta/callithrix_jacchus/dna/Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz -c 9483 -g http://ftp.ensembl.org/pub/release-104/gtf/callithrix_jacchus/Callithrix_jacchus.ASM275486v1.104.gtf.gz </code>
  - Users can update microbiome databases easily by one command line: **bash Update.sh -t [threads]**\
      For example:\
    <code>bash ~/MTD/Update.sh -t 20</code>
  - Users can modify the contaminant list (conta_ls.txt) in the MTD folder by adding the taxonomy ID of the microbe in the second column of the list and its name in the first column (optional).
  
# Licence
This software is freely available for academic users. Usage for commercial purposes is not allowed. Please refer to the LICENCE page.
