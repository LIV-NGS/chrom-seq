#!/usr/bin/python3

"""Pipeline for ATAC seq and ChIP seq

Usage:
 chrom-seq.py all-steps <fastq_dir> <home_dir> REF [-t <x>]
 chrom-seq.py align-pe <fqfile1> <fqfile2> REF [-t <x>] -o OUTPUT
 chrom-seq.py align-se <fqfile> REF [-t <x>] -o OUTPUT
 chrom-seq.py dedup-shift <bam_file> <home_dir> -t <x> [--atac|--chip]
 chrom-seq.py call-peaks <bam_file> <out_dir> <file_prefix> [--atac|--macs2|--epic2]
 chrom-seq.py (-h | --help)


Options:
  -h --help                 Show this screen.
   REF                      Full path to reference genome file.
  -t --threads=<x>          Number of threads to use for alignment [default: 4].
  -o --output=OUTPUT        Output file in BAM format for arguments align-pe, extension '.txt' when argument annotate (IMPORTANT: only use '_' in file names).
"""

from __future__ import print_function, division

import docopt
import errno
import os
import shutil
import sys
from subprocess import *

dir_list = ['picard.jar']


########################################################################################################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as ex:
        if ex.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def dir_names(homdir, *args):
    tdir = []
    ndir_list = dir_list
    try:
        for i in ndir_list:
            dc = Popen(["find", homdir, "-name", i], stdout=PIPE)
            work_dir = dc.stdout.read().rstrip()
            work_dir = work_dir.decode('utf-8')
            tdir.append(work_dir)
            dc.kill()
        tdir = tdir + list(args)
        # print(tdir)
        for i in tdir:
            if len(i) == 0:
                print("Cannot find %s, Please enter the correct parent directory for the mentioned tools: %s" % (
                    str(i), str(ndir_list[0])))
                sys.exit(-1)
        return tdir

    except CalledProcessError as cpe:
        print('Error %d.' % cpe.returncode)
        sys.exit(-1)


###########################
# ALIGN READS WITH Bowtie #
###########################


def align(fq1, fq2, ref, th, bamfile):
    bam_sort = bamfile.rstrip('.bam') + "_sorted.bam"
    align_mem = Popen(
        "bowtie2 --very-sensitive -X 2000  --no-mixed --no-discordant -p %d -x %s -1 %s -2 %s | samtools view -@ 5 -q 30 -bh -F 780 -f 2 - > %s" % (
            int(th), ref, fq1, fq2, bamfile),
        bufsize=-1, shell=True, executable='/bin/bash')
    align_mem.wait()
    samstat = Popen("samtools sort -@ 5 %s -o %s" % (bamfile, bam_sort), bufsize=-1, shell=True, executable='/bin/bash')
    samstat.wait()
    print("Indexing ...")
    Popen("samtools index %s" % (bam_sort), bufsize=-1, shell=True, executable='/bin/bash').wait()
    print("File is ready to deduplicate")


def alignse(fq1, ref, th, bamfile):
    bam_sort = bamfile.rstrip('.bam') + "_sorted.bam"
    align_mem = Popen(
        "bowtie2 --very-sensitive -p %d -x %s -U %s | samtools view -@ 5 -q 30 -bh -F 4 - > %s" % (
            int(th), ref, fq1, bamfile),
        bufsize=-1, shell=True, executable='/bin/bash')
    align_mem.wait()
    samstat = Popen("samtools sort -@ 5 %s -o %s" % (bamfile, bam_sort), bufsize=-1, shell=True, executable='/bin/bash')
    samstat.wait()
    print("Indexing ...")
    Popen("samtools index %s" % (bam_sort), bufsize=-1, shell=True, executable='/bin/bash').wait()
    print("File is ready to deduplicate")


########################################################################################################################


###################################
## FROM BAM FILE TO DEDUPLICATED ##
###################################


def dedup(bam, homdir, th):
    ndir = dir_names(homdir)

    bamfile = os.path.basename(bam)
    if not os.path.exists(''.join([bam.rstrip('.bam'), '_dedup.bam'])):
        print("Marking Duplicates...\n")
        dedp = Popen(
            "java -Xmx3g -jar %s MarkDuplicates I=%s O=%s M=%s REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate" \
            % (ndir[0], bamfile, ''.join([bamfile.rstrip('.bam'), '_dedup.bam']),
               ''.join([bamfile.rstrip('.bam'), '_dupstat.txt'])), bufsize=-1, shell=True, executable='/bin/bash')
        dedp.wait()
        print("Indexing ...")
        Popen("samtools index %s" % (''.join([bamfile.rstrip('.bam'), '_dedup.bam'])), bufsize=-1,
              shell=True, executable='/bin/bash').wait()

    if args['--atac']:
        if not os.path.exists(''.join([bam.rstrip('.bam'), '_shifted_final.bam'])):
            print("shifting reads...\n")
            ci = Popen(
                "java -Xmx3g -jar %s CollectInsertSizeMetrics I=%s O=%s H=%s " \
                % (ndir[0], ''.join([bamfile.rstrip('.bam'), '_dedup.bam']),
                   ''.join([bamfile.rstrip('.bam'), '_stat.txt']), ''.join([bamfile.rstrip('.bam'), '_hist.pdf'])),
                bufsize=-1, shell=True, executable='/bin/bash')
            ci.wait()

            dedp2 = Popen("alignmentSieve --numberOfProcessors %d  --ATACshift --bam %s --outFile %s" \
                          % (int(th), ''.join([bamfile.rstrip('.bam'), '_dedup.bam']),
                             ''.join([bamfile.rstrip('.bam'), '_dedup_shifted.bam'])), bufsize=-1, shell=True,
                          executable='/bin/bash')
            dedp2.wait()
            print("Sorting reads...\n")
            samstat2 = Popen("samtools sort -@ 5 %s -o %s" % \
                             (''.join([bamfile.rstrip('.bam'), '_dedup_shifted.bam']),
                              ''.join([bamfile.rstrip('.bam'), '_shifted_final.bam'])), bufsize=-1, shell=True,
                             executable='/bin/bash')
            samstat2.wait()
            print("Indexing ...")
            Popen("samtools index %s" % (''.join([bamfile.rstrip('.bam'), '_shifted_final.bam'])), bufsize=-1,
                  shell=True, executable='/bin/bash').wait()
            print("File is ready to Call peaks")
            mkdir_p("intermediate_files")
            intermed_files = [bam.rstrip('.bam') + '_dedup.bam', bam.rstrip('.bam') + '_dedup.bam.bai',
                              bam.rstrip('.bam') + '_dupstat.txt', bam.rstrip('.bam') + '_stat.txt',
                              ''.join([bamfile.rstrip('.bam'), '_dedup_shifted.bam'])]
            for intf in intermed_files:
                shutil.move(intf, "./intermediate_files")

    elif args['--chip']:
        print("File is ready to Call peaks")


def peaks(bam, out_dir, prename):
    if args['--atac']:
        print("Calling peaks ...\n")
        cp = Popen(
            "macs2 callpeak -t %s --format BAMPE --nomodel --call-summits --keep-dup all -B --SPMR -q 0.01 --outdir %s -n %s" % (
                bam, out_dir, prename), bufsize=-1, shell=True,
            executable='/bin/bash')
        cp.wait()
        print("\n\nGenerating bigwig tracks... \n")
        cp2 = Popen("bamCoverage --normalizeUsing BPM --binSize 10 --extendReads -p 35 -b %s -o %s" % (
            bam, ''.join([bam.rstrip('.bam'), '.bw'])), bufsize=-1, shell=True,
                    executable='/bin/bash')
        cp2.wait()
    elif args['--macs2']:
        if not os.path.exists(''.join([bam.rstrip('.bam'), '_input.bam'])):
            print("Input file is wrong or doesn't exists\n")
        else:
            print("Calling peaks with macs2 ...\n")
            cp = Popen(
                "macs2 callpeak -t %s -c %s --format BAMPE --keep-dup all -q 0.01 --outdir %s -n %s" % (
                    bam, bam.replace('.bam', '_input.bam'), out_dir, prename), bufsize=-1, shell=True,
                executable='/bin/bash')
            cp.wait()

        print("\n\nGenerating bigwig tracks... \n")
        cp2 = Popen("bamCoverage --normalizeUsing CPM --binSize 10 --extendReads -p 35 -b %s -o %s" % (
            bam, ''.join([bam.rstrip('.bam'), '.bw'])), bufsize=-1, shell=True,
                    executable='/bin/bash')
        cp2.wait()
    elif args['--epic2']:
        if not os.path.exists(''.join([bam.rstrip('.bam'), '_input.bam'])):
            print("Input file is wrong or doesn't exists\n")
        else:
            print("Calling peaks with epic2 ...\n")
            cp = Popen(
                "epic2 --treatment %s --control %s -fdr 0.05 --effective-genome-fraction 0.95 --guess-bampe --autodetect-chroms --output %s" % (
                    bam, bam.replace('.bam', '_input.bam'), prename), bufsize=-1, shell=True,
                executable='/bin/bash')
            cp.wait()
            mkdir_p("Peak_file")
            shutil.move(prename, "./Peak_file")

        print("\n\nGenerating bigwig tracks... \n")
        cp2 = Popen("bamCoverage --normalizeUsing CPM --binSize 10 --extendReads -p 35 -b %s -o %s" % (
            bam, ''.join([bam.rstrip('.bam'), '.bw'])), bufsize=-1, shell=True,
                    executable='/bin/bash')
        cp2.wait()


def all_steps(files_path, homdir, th, ref):
    try:
        bams_dir = [f for f in os.listdir(files_path) if
                    os.path.isfile(os.path.join(files_path, f)) and "R1_001.fastq.gz" in f]
    except OSError as ex:
        if ex.errno == errno.ENOENT:
            print("Enter correct directory path.")
            sys.exit(-1)
        else:
            raise

    bams = [f.replace("_R1_001.fastq.gz", ".bam") for f in bams_dir]
    bams_sorted = [f.replace("R1_001.fastq.gz", "sorted.bam") for f in bams_dir]
    bams_dedup = [f.replace("R1_001.fastq.gz", "sorted_dedup.bam") for f in bams_dir]
    bams_shifted = [f.replace("R1_001.fastq.gz", "dedup_shifted.bam") for f in bams_dir]
    bams_final = [f.replace("R1_001.fastq.gz", "shifted_final.bam") for f in bams_dir]

    for fq in bams_dir:
        align_mem = Popen(
            "bowtie2 --very-sensitive -X 2000  --no-mixed --no-discordant -p %d -x %s -1 %s -2 %s | samtools view -@ 5 -q 30 -bh -F 780 -f 2 - > %s" % (
                int(th), ref, fq, fq.replace("R1_001.fastq.gz", "R2_001.fastq.gz"), bams[bams_dir.index(fq)]),
            bufsize=-1, shell=True, executable='/bin/bash')
        align_mem.wait()
        samstat = Popen("samtools sort -@ 15 %s -o %s" % (bams[bams_dir.index(fq)], bams_sorted[bams_dir.index(fq)]),
                        bufsize=-1, shell=True,
                        executable='/bin/bash')
        samstat.wait()
        print("Indexing ...")
        Popen("samtools index %s" % (bams_sorted[bams_dir.index(fq)]), bufsize=-1, shell=True,
              executable='/bin/bash').wait()
        print("\n")
        print("#######################################")
        print("######### Removing duplicates.. #######")
        print("#######################################")
        print("\n")

        ndir = dir_names(homdir,dir_list)
        dedp = Popen(
            "java -Xmx3g -jar %s MarkDuplicates I=%s O=%s M=%s REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate" % (ndir[0], bams_sorted[bams_dir.index(fq)], bams_dedup[bams_dir.index(fq)],
               fq.replace("R1_001.fastq.gz", "_dupstat.txt")),
            bufsize=-1, shell=True, executable='/bin/bash')
        dedp.wait()
        print("Indexing ...")
        Popen("samtools index %s" % (bams_dedup[bams_dir.index(fq)]), bufsize=-1,
              shell=True, executable='/bin/bash').wait()

        print("\n")
        print("##################################")
        print("######### Shifting reads.. #######")
        print("##################################")
        print("\n")

        dedp2 = Popen("alignmentSieve --numberOfProcessors %d  --ATACshift --bam %s --outFile %s" % (int(th), bams_dedup[bams_dir.index(fq)], bams_shifted[bams_dir.index(fq)]),
                      bufsize=-1, shell=True, executable='/bin/bash')
        dedp2.wait()

        samstat2 = Popen("samtools sort -@ 5 %s -o %s" % \
                         (bams_shifted[bams_dir.index(fq)], bams_final[bams_dir.index(fq)]), bufsize=-1, shell=True,
                         executable='/bin/bash')
        samstat2.wait()
        print("Indexing ...")
        Popen("samtools index %s" % (bams_final[bams_dir.index(fq)]), bufsize=-1,
              shell=True, executable='/bin/bash').wait()

        ci = Popen(
            "java -Xmx3g -jar %s CollectInsertSizeMetrics I=%s O=%s H=%s " \
            % (ndir[0], bams_final[bams_dir.index(fq)],
               fq.replace("R1_001.fastq.gz", "_ISstat.txt"), fq.replace("_R1_001.fastq.gz", "IS.pdf")),
            bufsize=-1, shell=True, executable='/bin/bash')
        ci.wait()

        print("\n")
        print("#################################")
        print("######### Calling Peaks.. #######")
        print("#################################")
        print("\n")

        mkdir_p("intermediate_files")
        intermed_files = [bams_sorted[bams_dir.index(fq)], bams[bams_dir.index(fq)], bams_dedup[bams_dir.index(fq)],
                          bams_shifted[bams_dir.index(fq)], fq.replace("R1_001.fastq.gz", "_dupstat.txt"),
                          fq.replace("R1_001.fastq.gz", "_ISstat.txt")]
        for intf in intermed_files:
            shutil.move(intf, "./intermediate_files")

        cp = Popen(
            "macs2 callpeak -t %s --format BAMPE --nomodel --call-summits --keep-dup all -B --SPMR -q 0.01 --outdir %s -n %s" % (
                bams_final[bams_dir.index(fq)], fq.replace("_R1_001.fastq.gz", ""), fq.replace("R1_001.fastq.gz", "")),
            bufsize=-1, shell=True,
            executable='/bin/bash')
        cp.wait()

        print("\n\nGenerating bigwig tracks... \n")
        cp2 = Popen("bamCoverage --normalizeUsing BPM --binSize 10 --extendReads -p 35 -b %s -o %s" % (
            bams_final[bams_dir.index(fq)], fq.replace("_R1_001.fastq.gz", ".bw")), bufsize=-1, shell=True,
                    executable='/bin/bash')
        cp2.wait()

        print("\n")
        print("############################")
        print("######### FINISHED.. #######")
        print("############################")
        print("\n")


#####################################
#############    MAIN   #############
#####################################


if __name__ == '__main__':
    try:
        args = docopt.docopt(__doc__, version='chrom-seq-v1.0')

        if args['all-steps']:
            all_steps(args['<fastq_dir>'], args['<home_dir>'], args['--threads'], args['REF'])

        elif args['align-pe']:
            align(args['<fqfile1>'], args['<fqfile2>'], args['REF'], args['--threads'], args['--output'])

        elif args['align-se']:
            alignse(args['<fqfile>'], args['REF'], args['--threads'], args['--output'])

        elif args['dedup-shift']:
            dedup(args['<bam_file>'], args['<home_dir>'], args['--threads'])

        elif args['call-peaks']:
            peaks(args['<bam_file>'], args['<out_dir>'], args['<file_prefix>'])

    except docopt.DocoptExit:
        print("Operation Failed (Check Arguments) !! for help use: python chrom-seq.py -h or --help")
