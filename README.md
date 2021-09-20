# MTD: Meta-Transcriptome Detector
MTD is a software has two pipelines to detect and quantify microbiomes by analyzing bulk RNA-seq data and single-cell RNA-seq data, respectively. MTD is executed in Bash in GNU/Linux system. Users can easily install and run MTD using only one command and without requiring root privileges. The outputs (graphs, tables, count matrixes, etc.) are automatically generated and stored in the designated directory/folder defined by the user.
# Requirements
* 160 Gb RAM
* 550 Gb storage space
* Conda was installed
* GNU/Linux system with Bash
# Installation
1. Download directly or git clone MTD(git clone https://github.com/FEI38750/MTD.git) to the place you want to install. The full version of software (require 550Gb) will be installed in this MTD folder.
2. In termial, type **bash [path/to/MTD]/Install.sh -t [threads] -p [path/to/conda]**\
For example:
<pre><code>bash ~/MTD/Install.sh -t 20 -p ~/miniconda3
</code></pre>
## Notes
* Installation may take 1-2 days.
* If conda hasn't been installed in your system, please use the code below to install the conda:
<addr><pre><code>URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
curl $URL > miniconda-installer.sh
bash miniconda-installer.sh -b</code></pre>
Conda then will be installed in your home directory, such as path: ~/miniconda3
* Tips: file management software such as FileZilla (https://filezilla-project.org/download.php?show_all=1) can help you to manage your files on HPC/server.
# Run MTD
## Bulk RNA-seq
1. Store all the paired-end fastq files (accepted: fastq, fastq.gz, fq, fq.gz) to be analyzed in a folder, subfolders for each sample are accepted.
  The paired fastq files must be named starting with the sample name followed by "_1" and "_2". For example, sample1_1.fq.gz and sample1_2.fq.gz are paired-end fastq files for sample1.
2. Prepare the samplesheet.csv. You can copy and modify the one in MTD folder.
  ![image1](https://github.com/FEI38750/MTD/blob/main/Img/Tutorial1.jpg)
3. Put samplesheet.csv in the same folder as the fastq files.
  ![image2](https://github.com/FEI38750/MTD/blob/main/Img/Tutorial2.jpg)
3. In termial, type **bash [path/to/MTD]/MTD.sh -i [path/to/samplesheet.csv] -o [path/to/output_folder] -h [host species taxonomy ID] -t [threads]**\
Host species taxonomy ID: human:9606, mouse:10090, rhesus monkey:9544\
For example:
<pre><code>bash ~/MTD/MTD.sh -i ~/raw_data/samplesheet.csv -o ~/MTD_output/ -h 9544 -t 20</code></pre>
## Single-cell RNA-seq
In termial, type **bash [path/to/MTD]/MTD_singleCell.sh -i [path/to/Input_folder] -o [path/to/Output_folder] -h [Host species taxonomy ID] -t [Threads] -p [Platform] -d [prime Direction] -c [path/to/Cell_barcode_file.whitelist.txt]**\
  Single cell RNAseq platform(-p): enter 1 for 10x or 2 for Dropseq platform\
  For example:
<pre><code>bash ~/MTD/MTD_singleCell.sh -i ~/scRNAseq_rawData/ -o ~/output/ -h 10090 -t 20 -p 1 -d 3 -c ~/scRNAseq_rawData/SRR12345678.whitelist.txt</code></pre>
### Notes
* 10x and Dropseq use paired end sequence. The second fastq file contains transcript's sequences (e.g., 98bp length). Please run MTD_singleCell.sh for the second reads.
* -i [path/to/Input_folder] contains the single-cell RNAseq raw data in fastq format (accepted: fastq, fastq.gz, fq, fq.gz). For example: ~/scRNAseq_rawData/ contains Scsample1_1.fastq and Scsample1_2.fastq for sample Scsample1.
* The paired fastq files must be named starting with the sample name followed by "_1" and "_2". For example, Scsample1_1.fastq and Scsample1_2.fastq are paired end fastq files for Scsample1. Scsample1_1.fastq contains barcodes (e.g., 26bp length) and Scsample1_1.fastq contains transcript's sequences (e.g., 98bp length).
* The whitelist.txt contains the cell barcodes gotten from the normal single-cell RNAseq analysis for the host reads. Each cell barcode occupies one line.
  ![whitelist_example](https://github.com/FEI38750/MTD/blob/main/Img/whitelist_img.png)
# Licence
This software is freely available for academic users. Usage for commercial purposes is not allowed. Please refer to the LICENCE page.
