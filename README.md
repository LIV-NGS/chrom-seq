# chrom-seq**
Python pipeline for ATAC-seq and ChIP-seq analysis
**
Install following tools and setup in your path (except picard.jar):

bowtie2
MACS2
PICARD
samtools


**usage:**

chrom-seq.py -h 

**example:
**
chrom-seq.py align-pe sample_read1.fastq.gz sample_read2.fastq.gz ~/Desktop/ref/Bowtie2Index -t 3 -o test.bam

chrom-seq.py dedup-shift test_sorted.bam ~/Desktop/tools/ -t 6 --atac

~/Desktop/work_hpi/python/chrom-seq.py call-peaks test_sorted_dedup.bam test peaks --atac
