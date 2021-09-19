# MTD: Meta-Transcriptome Detector
The MTD has two pipelines to detect and quantify microbiomes by analyzing bulk RNA-seq data and single-cell RNA-seq data, respectively. MTD is executed in GNU/Linux system. Users can easily install and run MTD using only one command and without requiring root privileges. The outputs (graphs, tables, count matrixes, etc.) are automatically generated and stored in the designated directory/folder defined by the user.
# Installation
Requirements:
* 160 Gb RAM
* 550 Gb storage space
* Conda was installed
# Install
1. Download directly or git clone MTD(git clone https://github.com/FEI38750/MTD.git) to the place you want to install. The full version of software (require 550Gb) will be installed in this MTD folder.
2. In termial, type **bash [path/to/MTD]/Install.sh -t [threads] -p [path/to/conda]**\
For example:
<pre><code>bash ~/MTD/Install.sh -t 20 -p ~/miniconda3
</code></pre>
Note: Installation may take 1-2 days.
