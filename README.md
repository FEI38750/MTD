# MTD: Meta-Transcriptome Detector
The MTD has two pipelines to detect and quantify microbiomes by analyzing bulk RNA-seq data and single-cell RNA-seq data, respectively. MTD is executed in GNU/Linux system. Users can easily install and run MTD using only one command and without requiring root privileges. The outputs (graphs, tables, count matrixes, etc.) are automatically generated and stored in the designated directory/folder defined by the user.
# Installation
Requirements:
* 160 Gb RAM
* 550 Gb storage space
# Install
1. Download directly or git clone MTD(git clone https://github.com/FEI38750/MTD.git)
2. Open the setup.sh by a text editor. Modify parameters: MTD_path=(path of the MTD folder), threads=(number of threads you want to use for installation), conda_path=(path of miniconda/anaconda on your system). Then save and run setup.sh.
<pre><code>bash setup.sh
</code></pre>

