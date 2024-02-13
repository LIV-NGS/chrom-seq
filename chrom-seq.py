#!/usr/bin/python3

"""Pipeline for ATAC seq and ChIP seq

Usage:
 chrom-seq.py align-pe <fqfile1> <fqfile2> REF [-t <x>] -o OUTPUT
 chrom-seq.py dedup-shift <bam_file> <home_dir> -t <x> [--atac|--chip]
 chrom-seq.py call-peaks <bam_file> <out_dir> <dir_name> [--atac|--chip]
 chrom-seq.py (-h | --help)


Options:
  -h --help                 Show this screen.
   REF                      Full path to reference genome file.
  -t --threads=<x>          Number of threads to use for alignment [default: 2].
  -o --output=OUTPUT        Output file in BAM format for arguments align-pe, extension '.txt' when argument annotate (IMPORTANT: only use '_' in file names).
"""


from __future__ import print_function, division
from  subprocess import *
#from statistics import *
#from numpy import *
import sys, os, shutil, csv, time, docopt, errno, itertools, gzip, re




#pic_cmd = ["AddOrReplaceReadGroups", "SortSam", "MarkDuplicates"]
#pic_arg = ["I=", "O=", "VALIDATION_STRINGENCY=LENIENT", "M=", 'SORT_ORDER=', 'CREATE_INDEX=']
#dir_list = ['picard.jar']

r = re.compile('(\d+)')

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
        #print(tdir)
        for i in tdir:
            if len(i) == 0:
                print("Cannot find %s, Please enter the correct parent directory for the mentioned tools: %s" % (
                str(i), str(ndir_list[0])))
                sys.exit(-1)
        return tdir

    except CalledProcessError as cpe:
        print('Error %d.' % cpe.returncode)
        sys.exit(-1)



########################
# ALIGN READS WITH BWA #
########################


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



########################################################################################################################


#############################################################
## FROM BAM FILE TO DEDUPLICATED AND RECALIBRATED BAM FILE ##
#############################################################


def add_RG(bam, homdir,th):
    ndir = dir_names(homdir)

    bamfile = os.path.basename(bam)
    if not os.path.exists(''.join([bam.rstrip('.bam'), '_dedup.bam'])):
        print("Marking Duplicates...\n")
        dedp = Popen("java -Xmx3g -jar %s MarkDuplicates I=%s O=%s M=%s REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate" \
                     % (ndir[0],bamfile, ''.join([bamfile.rstrip('.bam'), '_dedup.bam']), ''.join([bamfile.rstrip('.bam'), '_dupstat.txt'])), bufsize=-1,shell=True, executable='/bin/bash')
        dedp.wait()
        print("Indexing ...")
        Popen("samtools index %s" % (''.join([bamfile.rstrip('.bam'), '_dedup.bam'])), bufsize=-1,
              shell=True, executable='/bin/bash').wait()

    if args['--atac']:
        if not os.path.exists(''.join([bam.rstrip('.bam'), '_shifted_final.bam'])):
            print("shifting reads...\n")
            ci = Popen(
                "java -Xmx3g -jar %s CollectInsertSizeMetrics I=%s O=%s H=%s " \
                % (ndir[0], ''.join([bamfile.rstrip('.bam'), '_dedup.bam']),''.join([bamfile.rstrip('.bam'), '_stat.txt']) ,''.join([bamfile.rstrip('.bam'), '_hist.pdf'])),
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
    elif args['--chip']:
        print("File is ready to Call peaks")




def peaks(bam,out_dir,prename):
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
    elif args['--chip']:
        if not os.path.exists(''.join([bam.rstrip('.bam'), '_input.bam'])):
            print("Input file is wrong or doesn't exists\n")
        else:
            print("Calling peaks ...\n")
            cp = Popen(
                "macs2 callpeak -t %s -c %s --format BAMPE --keep-dup all --broad  -q 0.05 --broad-cutoff 0.05 --outdir %s -n %s" % (
                    bam, bam.replace('_dedup.bam', '_input.bam'), out_dir, prename), bufsize=-1, shell=True,
                executable='/bin/bash')
            cp.wait()

        print("\n\nGenerating bigwig tracks... \n")
        cp2 = Popen("bamCoverage --normalizeUsing CPM --binSize 10 --extendReads -p 35 -b %s -o %s" % (
            bam, ''.join([bam.rstrip('.bam'), '.bw'])), bufsize=-1, shell=True,
                    executable='/bin/bash')
        cp2.wait()






#####################################
#############    MAIN   #############
#####################################


if __name__ == '__main__':
    try:
        args = docopt.docopt(__doc__, version='chrom-seq-v1.0')

        if args['align-pe']:
            align(args['<fqfile1>'], args['<fqfile2>'], args['REF'], args['--threads'], args['--output'])


        elif args['dedup-shift']:
            add_RG(args['<bam_file>'], args['<home_dir>'], args['--threads'])

        elif args['call-peaks']:
            peaks(args['<bam_file>'], args['<out_dir>'], args['<dir_name>'])

    except docopt.DocoptExit:
        print("Operation Failed (Check Arguments) !! for help use: python ATAC-seq.py -h or --help")
