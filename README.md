[![R](https://ziadoua.github.io/m3-Markdown-Badges/badges/R/r1.svg)](https://www.r-project.org/)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/LIV-NGS/RNA-seq-Report/graphs/commit-activity)
# chrom-seq


**Python pipeline for ATAC-seq and ChIP-seq analysis**

Install following tools and setup in your path (except picard.jar):

bowtie2: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.3/

MACS2: https://pypi.org/project/MACS2/

PICARD: https://broadinstitute.github.io/picard/

samtools: http://www.htslib.org/download/ 

docopt: https://github.com/docopt/docopt


**usage:**

chrom-seq.py -h 

**example:**

```
chrom-seq.py all-steps /home/user/fastq_dir /home ~/Desktop/ref/Bowtie2Index -t 34

chrom-seq.py align-pe sample_read1.fastq.gz sample_read2.fastq.gz ~/Desktop/ref/Bowtie2Index -t 3 -o test.bam

chrom-seq.py dedup-shift test_sorted.bam ~/Desktop/tools/ -t 6 --atac

chrom-seq.py call-peaks test_sorted_dedup.bam test peaks --atac
```
________________________________________________________________________________________________

**Differential ATAC-seq report**
R markdown document to creat html reports for Differential ATAC-seq.

Following are the steps to use the R markdown script.

Install following tools and setup in your path (except picard.jar):

deeptools: https://deeptools.readthedocs.io/en/latest/

HOMER (annotatePeaks.pl): http://homer.ucsd.edu/homer/ngs/annotation.html

Type echo $PATH in your terminal and initialize path at the start of the script (install all the packages too):

1. https://github.com/ColeWunderlich/soGGi


```
### Replace "..."  with output of echo $PATH
Sys.setenv(PATH = paste(old_path, "...", sep = ":")) 
###Run
rmarkdown::render("PATH-to-the-rmd-script",output_dir = "full_path_to_output_dir")
```
