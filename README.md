# MTD: Meta-Transcriptome Detector
MTD is a software has two pipelines to detect and quantify microbiomes by analyzing bulk RNA-seq data and single-cell RNA-seq data, respectively. MTD is executed in Bash in GNU/Linux system. Users can easily install and run MTD using only one command and without requiring root privileges. The outputs (graphs, tables, count matrixes, etc.) are automatically generated and stored in the designated directory/folder defined by the user.
# Requirements
* 160 Gb RAM
* 550 Gb storage space
* Conda was installed
# Installation
1. Download directly or git clone MTD(git clone https://github.com/FEI38750/MTD.git) to the place you want to install. The full version of software (require 550Gb) will be installed in this MTD folder.
2. In termial, type **bash [path/to/MTD]/Install.sh -t [threads] -p [path/to/conda]**\
For example:
<pre><code>bash ~/MTD/Install.sh -t 20 -p ~/miniconda3
</code></pre>
# Notes
* Make the script executable with command chmod +x <fileName> before installation, when nessary. For example: 
<addr><pre><code>chmod +x ~/MTD/Install.sh</code></pre>
* Installation may take 1-2 days.
* If conda hasn't been installed in your system, please use the code below to install the conda:
<addr><pre><code>URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
curl $URL > miniconda-installer.sh
bash miniconda-installer.sh -b</code></pre>
Conda will be installed in your home directory, such as path: ~/miniconda3\
* File management software such as FileZilla (https://filezilla-project.org/download.php?show_all=1) can help you to manage your files on HPC.

