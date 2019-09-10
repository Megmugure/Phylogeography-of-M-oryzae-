#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import csv
import argparse
import logging
import subprocess
import logging.config
import getpass
import time
from datetime import datetime

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
logfile = os.path.splitext(os.path.basename(os.path.abspath(__file__)))[0]+'.log'
logging.basicConfig(filename=logfile, filemode='w', level=logging.INFO,
format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.getLogger().addHandler(logging.StreamHandler())

# set start time format
format = "%a %b %d %H:%M:%S %Y"
start_time = datetime.now()
s = start_time.strftime(format)
logging.info("\nStarted analysis on: " +s+  "\n")


# this script performs alignment based SNP calling
# reference genome => ftp://ftp.ensemblgenomes.org/pub/release-43/fungi/fasta/magnaporthe_oryzae/dna/Magnaporthe_oryzae.MG8.dna.toplevel.fa.gz


def run_fastqc(indir, outdir, t=0):
    """
    perform quality control on NGS data using FastQC (v0.11.7)
    parameters:
        - indir (str): directory having raw sqeuencing data (.fq, .fastq, .fastq.gz, .fq.gz)
        - outdir (str): output directory for the FastQC analysis
        - t (str): number of threads/cores
    returns:
        - reads_dict (str): hash table/dictioanry with isolate and reads

    """

    fastq_ext = 'fastq fq gz'.split()
    isolates_dict = dict()
    fastq_dict = dict()
    for root, dirs, files in os.walk(indir):
        for dir in dirs:
            if dir.startswith('Isolate'):
                isolate = dir
                path = os.path.join(root, dir)
                if not isolate in isolates_dict:
                    isolates_dict[isolate] = path
    for isolate, path in isolates_dict.items():
        for file in sorted(os.listdir(path)):
            f = os.path.splitext(file)[1][1:]
            if not f in fastq_ext:
                continue
            else:
                fq = os.path.join(path, file)
                if not isolate in fastq_dict:
                    fastq_dict[isolate] = [fq]
                else:
                    fastq_dict[isolate].append(fq)


    for isolate, reads in fastq_dict.items():
        for read in reads:
            base = os.path.basename(read)
            html = os.path.join(outdir, base.rsplit('.')[0]+'_fastqc.html')
            zip = os.path.join(outdir, base.rsplit('.')[0]+'_fastqc.zip')
            if os.path.exists(html) and os.path.exists(zip):
                continue
            else:
                cmd = 'fastqc --nogroup {} -o {} --threads {}'.format(read, outdir, t)
                try:
                    logging.info("performing QC on: {}\t{} threads".format(os.path.basename(read), t))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception:
                    logging.error("\nError occurred when performing QC.\n")
                    logging.info("\ncommand that was running {}\n".format(cmd))

    return fastq_dict


def build_genome_index(genome_fa):
    """
    Build genome index using BWA v0.7.15

    parameters:
        - genome_fa (str) : reference file in FASTA format
        - genome_index (str): genome index
    returns:
        - index_base (str): base index for the genome
    """

    index_base = os.path.splitext(os.path.abspath(genome_fa))[0]
    # -p prefix of the index [same as fasta name]
    cmd = 'bwa index -p {} {}'.format(index_base, os.path.abspath(genome_fa))
    idx = os.path.splitext(os.path.abspath(genome_fa))[0]+'.bwt'
    if os.path.exists(idx):
        logging.info("\n{} reference index exists\n".format(os.path.basename(genome_fa)))
    else:
        try:
            logging.info("\n{} reference index does not exist, building index\n\
            ".format(os.path.basename(genome_fa)))
            p = subprocess.check_call(cmd, shell=True)
        except Exception as e:
            logging.error("Error occured when building genome index for {}".format(genome_fa))
            logging.error("Error: {}".format(e))
    return index_base



def align_reads(reads_dict, index_base, outdir, t=0):
    """
    map/align reads to the reference genome using BWA  v0.7.15
    use BWA so that read group information can be added during
    the alignment stage, without requiring a separate step with
    Picard's AddOrReplaceReadGroups tool.

    parameters:
        - reads_dict (str): hash table/dictionary with isolate and reads
        - outdir (str): directory to store the output of the mapping
        - index_base (str) : base index for the genome
        - t (int): number of threads/cores
    returns:
        - outdir (str): directory having alignments

    """

    for isolate, reads in reads_dict.items():
        if len(reads) == 2:
            flowcell = reads[0].rsplit('_')[2]
            lane = reads[0].rsplit('_')[1].strip('L')
            sample = isolate
            read_group = repr('@RG\tID:{}.{}\tSM:{}'.format(flowcell, lane, sample))
            outname = isolate+'_'+reads[0].rsplit('_')[2]+'_'+reads[0].rsplit('_')[-1]
            pe_name = outname.rsplit('_', 1)[0]

            sam_out = pe_name+''+'.sam'
            samfilename = os.path.join(outdir, sam_out)
            r1, r2 = reads[0], reads[1]

            # -R (str) read group header line such as '@RG\tID:{FLOWCELL}.{LANE}\tPU:{FLOWCELL_BARCODE}.{LANE}.{SAMPLE}\tSM:{SAMPLE}\tPL:{PLATFORM}\tLB{LIBRARY}
            cmd = "bwa mem -t {} -R {} {} {} {} > {}".format(t, read_group, index_base, r1, r2, samfilename)
            if os.path.exists(samfilename) and os.path.getsize(samfilename) > 0:
                print("{} & {} already mapped, check {}".format(
                os.path.basename(r1), os.path.basename(r2), samfilename))
                continue
            else:

                try:
                    logging.info("mapping: {}\t{}\t{}\t{}".format(os.path.basename(r1), os.path.basename(r2), os.path.basename(index_base), t))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occured when mapping {} & {} to {}".format(os.path.basename(r1), os.path.basename(r2), os.path.basename(index_base)))
                    logging.error("Error: {}".format(e))
    return outdir

def sortIndexsam(indir):
    """
    Sort and index the alignment and store as a binary file (BAM)
    sort SAM file by coordinate position (required for downstream analyses)
    and output the file in BAM format with PicardTools SortSam command

    parameters:
        - indir (str): directory having the alignments

    returns:
        - indir (str): directory having the alignments
    """

    tool = 'SortSam'

    for f in sorted(os.listdir(indir)):
        if f.endswith('.sam') or f.endswith('.bam'):
            base = f.rsplit('.')[0]
            samfile = os.path.join(indir, f)
            bamfile = os.path.join(indir, base.rsplit('.')[0])+'.bam'

            sorted_bam = os.path.splitext(bamfile)[0]+'.sorted.bam'
            sorted_bam = os.path.join(indir, sorted_bam)

            if os.path.exists(sorted_bam):
                continue
            else:
                # SORT_INDEX Sort order of output file  Required. Possible values: {unsorted, queryname, coordinate, duplicate}
                sort_order = 'coordinate'
                create_index = 'true'
                cmd = 'java -jar /export/apps/picard/2.8.2/picard.jar {} I={} \
                O={} SORT_ORDER={} \
                CREATE_INDEX={}'.format(tool, samfile, sorted_bam, sort_order,
                create_index)
                try:
                    logging.info("Sorting SAM file {} to {}".format(os.path.basename(samfile), os.path.basename(sorted_bam)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occured when sorting {}".format(os.path.basename(samfile)))
                    logging.error("Error: {}".format(e))
    return indir

def collect_alignmentMetrics(indir, reference):
    """
    calculate the alignment metrics for samples using Picard tool (v2.8.2)
    CollectAlignmentSummaryMetrics

    parameters:
        - indir (str): directory having the alignments
        - reference (str): reference genome sequence in FASTA format
    returns:
        - indir (str): directory having the alignments
    """

    tool = 'CollectAlignmentSummaryMetrics'

    for f in sorted(os.listdir(indir)):
        if f.endswith('.bam'):
            base = f.rsplit('.')[0]
            bamfile = os.path.join(indir, f)
            metrics_file = os.path.join(indir, base.rsplit('.')[0])+'.metrics.txt'

            if os.path.exists(metrics_file):
                print("alignment metrics already collected for: {}".format(os.path.basename(bamfile)))
                continue
            else:
                cmd = "java -jar /export/apps/picard/2.8.2/picard.jar {}\
                I={} R={} METRIC_ACCUMULATION_LEVEL=SAMPLE \
                METRIC_ACCUMULATION_LEVEL=READ_GROUP \
                O={}".format(tool, bamfile, reference, metrics_file)
                try:
                    logging.info("collecting alignment metrics for {}".format(os.path.basename(bamfile)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occured when collecting metrics for {}".format(os.path.basename(bamfile)))
                    logging.error("Error: {}".format(e))
    return indir


def deduplicate_bam(indir, outdir, t):
    """
    identify any duplicate sequences from the same DNA fragment BAM files
    This occur due to sample preparation (e.g. during PCR) or incorrect optical
    cluster identification during sequencing

    using Picard tool (v2.8.2) MarkDuplicates

    You can also merge all BAM files originating frrom the same sample at this
    point. (see options)

    parameters:
        - indir (str): directory having filtered BAM alignments
        - outdir (str): directory to store the deduplicated BAM alignments
        - t (int): number of threads/cores

    returns:
        - outdir (str): directory to store deduplicated BAM alignments
    """

    '''program options'''
    tool = 'MarkDuplicates'
    # it is not recommended to actually remove the duplicate sequences from the
    # file, but simply to mark the flags appropriately in the BAM file,
    # so that those sequences are ignored downstream.
    # If using tools other than those used in this pipleine, make sure they can
    # identify these flags.
    remove_dup = 'false'
    tagging_policy = 'All'

    for f in sorted(os.listdir(indir)):
        if f.endswith('.bam'):
            base = f.rsplit('.')[0]
            infile = os.path.join(indir, f)
            outfile = os.path.join(outdir, base)+'_dedup.bam'
            metrics = os.path.join(outdir, base)+'_dedup_metrics.txt'

            if os.path.exists(outfile):
                print("{} already exists!".format(outfile))
                continue
            else:
                cmd = "java -jar /export/apps/picard/2.8.2/picard.jar {}\
                INPUT={} OUTPUT={} METRICS_FILE={} REMOVE_DUPLICATES={} \
                TAGGING_POLICY={}".format(tool, infile, outfile, metrics,
                remove_dup, tagging_policy)
                try:
                    logging.info("identifying duplicates in {}".format(infile))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when removing duplicate: {}".format(e))
    return outdir

def validate_bamfile(indir):
    """
    validate BAM file to make sure there were not issues or
    mistakes associated with previous analyses.
    using Picard ValidateSamFile tool (v2.8.2)
    look at the errors you may want to ignore , in this case, I choose to
    ignore the MISSING_PLATFORM_VALUE error. I presume platform was Illumina but
    since I did not include the information while aligning/mapping the reads
    in the alignment step using the -R (READ GROUP) option, I am bound to
    experience such an error.
    Refer to: https://software.broadinstitute.org/gatk/documentation/article.php?id=7571

    ///This step is not mandatory but is important for downstream processes///

    parameters:
        - indir (str): directory having deduplicated BAM files

    returns:
        - indir (str): directory having deduplicated BAM files
    """
    tool = 'ValidateSamFile'

    for fn in sorted(os.listdir(indir)):
        if fn.endswith('.sorted.bam'):
            isolate = fn.rsplit('_')[0]
            validate_file = os.path.join(indir, isolate)+'.validate.txt'
            infile = os.path.join(indir, fn)

            if os.path.exists(validate_file):
                continue
            else:
                # ignore errors
                # IGNORE=type
                error_1 = 'MISSING_PLATFORM_VALUE'
                cmd = 'java -jar /export/apps/picard/2.8.2/picard.jar {} I={}\
                 O={} MODE=SUMMARY IGNORE={}\
                '.format(tool, infile, validate_file, error_1)

                try:
                    logging.info("validating the BAM file {}".format(infile))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when validating the file: {}".format(os.path.basename(infile)))
                    logging.error("Error that occurred: {}".format(e))
    return indir


def create_reference_index(reference):
    """
    GATK uses two files to access and safety check access to the reference files:
    a .dict dictionary of the contig names and sizes and a .fai fasta index file
    to allow efficient random access to the reference bases.
    You have to generate these files in order to be able to use a Fasta file as
    reference.

    uses:
        - CreateSequenceDictionary.jar from Picard for dict
        - samtools (v1.8) faidx


    parameters:
        - reference (str): reference FASTA format file

    returns:
        - dict (str): dictionary of the contig names and sizes
        - index (str): fasta index file
    """
    index_base = os.path.splitext(os.path.abspath(reference))[0]
    dict = index_base+'.dict'
    fai = reference+'.fai'
    tool = 'CreateSequenceDictionary'

    if os.path.exists(dict):
        print("{} dictionary already exists".format(os.path.basename(dict)))

    else:
        cmd = 'java -jar /export/apps/picard/2.8.2/picard.jar {} R={} O={}\
        '.format(tool, reference, dict)
        try:
            logging.info("creating a dictionary for the reference {}".format(os.path.basename(reference)))
            p = subprocess.check_call(cmd, shell=True)
        except Exception as e:
            logging.error("Error occurred when creating a dict for the reference {}".format(os.path.basename(reference)))
            logging.error("Error that occurred: {}".format(e))


    if os.path.exists(fai):
        print("{} index already exists".format(os.path.basename(fai)))

    else:
        cmd = 'samtools faidx {}'.format(reference)
        try:
            logging.info("creating an index for the reference {}".format(os.path.basename(reference)))
            p = subprocess.check_call(cmd, shell=True)
        except Exception as e:
            logging.error("Error occurred when creating index for the reference {}".format(os.path.basename(reference)))
            logging.error("Error that occurred: {}".format(e))
    return reference

def target_indel_intervals(indir, reference, t):
    """
    define intervals to target for local realignment using RealignerTargetCreator
    tool in GATK (v3.7.0)

    parameters:
        - indir (str): directory having deduplicated sorted alignments in BAM format
        - reference (str): reference genome sequence in FASTA format
        - outdir (str): directory to store the output (intervals files)
        - t (int): number of threads

    returns:
        - outdir (str): directory to store the output (intervals files)
    """
    # define intervals to target for local realignment
    # java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 4
    #-R ~/genomes/ensembl/Magnaporthe_oryzae/assembly/Magnaporthe_oryzae.MG8.dna.toplevel.fa
    #-I /var/scratch/jjuma/testdata_output/picard_output/Isolate71_wHAXPI092787-141_dedup.sorted.bam
    # -o /var/scratch/jjuma/testdata_output/picard_output/Isolate71_wHAXPI092787-141_realign.intervals

    tool = 'RealignerTargetCreator'

    for f in sorted(os.listdir(indir)):
        if f.endswith('.sorted.bam'):
            base = f.rsplit('_', 1)[0]
            infile = os.path.join(indir, f)
            outfile = os.path.join(indir, base)+'_realign.intervals'

            if os.path.exists(outfile):
                continue
            else:
                cmd = 'java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T\
                {} -nt {} -R {} -I {} -o {}'.format(tool, int(t/t),
                reference, infile,outfile)
                try:
                    logging.info("Defining intervals to target for local realignment in {}".format(os.path.basename(infile)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when targeting intervals in {}".format(os.path.basename(infile)))
                    logging.error("Error that occurred: {}".format(e))
    return indir

def realign_indels(indir, reference, t):
    """
    perform local realignment of reads around indels using IndelRealigner
    GATK (v3.7.0)

    parameters:
        - indir (str): directory having deduplicated sorted alignments in BAM format and indels intervals per sample
        - reference (str): reference genome sequence in FASTA format
        - t (int): number of threads

    returns:
        - indels_dir (str): directory having samples/isolates' realignment intervals
    """

    # perform local realignment of reads around indels
    # java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner
    #-R ~/genomes/ensembl/Magnaporthe_oryzae/assembly/Magnaporthe_oryzae.MG8.dna.toplevel.fa
    #-targetIntervals /var/scratch/jjuma/testdata_output/picard_output/Isolate71_wHAXPI092787-141_realign.intervals
    #-I /var/scratch/jjuma/testdata_output/picard_output/Isolate71_wHAXPI092787-141_dedup.sorted.bam
    #-o /var/scratch/jjuma/testdata_output/picard_output/Isolate71_wHAXPI092787-141_dedup.sorted.realigned.bam

    tool = 'IndelRealigner'
    maxreads = 2147483647 # maximum reads allowed at an interval for realignment

    isolates_dict = dict()
    for f in sorted(os.listdir(indir)):
        if not f.endswith('.sorted.bam') and not f.endswith('.intervals'):
            continue
        else:
            if f.endswith('.sorted.bam'):
                base = f.rsplit('_', 1)[0]
                infile = os.path.join(indir, f)
                isolates_dict[base] = [infile]
            else:
                indel_file = os.path.join(indir, f)
                if base in isolates_dict:
                    isolates_dict[base].append(indel_file)



    for isolate, fns in isolates_dict.items():
        outfile = os.path.splitext(fns[0])[0]+'.realigned.bam'
        if os.path.exists(outfile):
            continue
        else:
            cmd = 'java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T\
            {} -nt {} -R {} -targetIntervals {} -maxReads {} -I {} -o {}\
            '.format(tool, int(t/t), reference, fns[1], maxreads, fns[0], outfile)
            try:
                logging.info("Realigning {} around {} intervals".format(os.path.basename(fns[0]), os.path.basename(fns[1])))
                p = subprocess.check_call(cmd, shell=True)
            except Exception as e:
                logging.error("Error occurred when realigning {} using {} intervals".format(os.path.basename(fns[0]), os.path.basename(fns[1])))
                logging.error("Error that occurred: {}".format(e))
    return indir

def fix_matepair_info(indir):
    """
    Fix mate pair information in BAM using FixMateInformation tool in Picard
    (v2.8.2)

    parameters:
        - indir (str): directory having realigned BAM alignments

    returns:
        - indir (str): directory having sorted realigned BAM alignments
    """
    tool = 'FixMateInformation'
    so = 'coordinate'
    create_index = 'true'

    for f in sorted(os.listdir(indir)):
        if not f.endswith('.realigned.bam'):
            continue
        else:
            infile = os.path.join(indir, f)
            outfile = os.path.join(indir, os.path.splitext(f)[0])+'.fixmate.bam'
            if os.path.exists(outfile):
                continue
            else:
                cmd = 'java -jar /export/apps/picard/2.8.2/picard.jar {} INPUT={}\
                OUTPUT={} SO={} CREATE_INDEX={}'.format(tool, infile, outfile, so, create_index)
                try:
                    logging.info("Fixing mate pair information in {} ".format(os.path.basename(infile)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when fixing mate pait info in {}".format(os.path.basename(infile)))
                    logging.error("Error that occurred: {}".format(e))
    return indir


def call_variants_UnifiedGenotyper(indir, reference, outdir, t):
    """
    call variants using GATK (v3.7.0) UnifiedGenotyper - SNP calling site-by-site

    parameters:
        - indir (str): directory having sorted deuplicated alignments in binary format
        - reference (str): reference genome sequence in FASTA format
        - outdir (str): directory to store the outputs
        - t (int): number of threads (cores)

    returns:
        - outdir (str): directory to store the outputs
    """

    tool = 'UnifiedGenotyper'
    variant_qual = 30

    for f in sorted(os.listdir(indir)):
        if f.endswith('.fixmate.bam'):
            base = f.rsplit('.')[0]
            infile = os.path.join(indir, f)
            vcf_file = os.path.join(outdir, base.rsplit('_')[0])+'.vcf'



            if os.path.exists(vcf_file):
                logging.info("variant calling on {} already done, check! {}\
                ".format(os.path.basename(infile), os.path.basename(vcf_file)))

            else:
                cmd = "java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T \
                {} -R {} -I {} -stand_call_conf {}\
                -o {}".format(tool, reference, infile, variant_qual, vcf_file)
                try:
                    logging.info("performing variant calling on {} \
                    ".format(os.path.basename(infile)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when performing variant calling")
                    logging.error("Error that occurred: {}".format(e))
    return outdir


def call_variants_HaplotypeCaller(indir, reference, outdir):
    """
    Call variants using  HaplotypeCaller program from the GATK (v3.7.0) pipeline.
    Calling variants with HaplotypeCaller is essentially a two-step process
    (similar to indel realignment). First, you call genotypes individually for
    each sample. Second, you perform joint genotyping across samples to produce
    a multi-sample VCF call-set. The advantage to this strategy is that the most
    computationally intensive step, calling genotypes for each sample,
    only needs to be performed once, even if additional samples will be added
    later. The joint genotyping, which is less computationally intensive, can
    be performed as many times as needed as individuals may be added to the
    dataset.

    The output for this program will be a GVCF, which has raw, unfiltered SNP
    and indel calls for all sites, variant or invariant, unlike a typical VCF
    file.

    This is specified by the --emitRefConfidence GVCF

    parameters:
        - indir (str): directory having sorted deuplicated alignments in binary format
        - reference (str): reference genome sequence in FASTA format
        - outdir (str): directory to store the outputs

    returns:
        - outdir (str): directory to store the outputs
    """

    for f in sorted(os.listdir(indir)):
        if f.endswith('.fixmate.bam'):
            base = f.rsplit('.')[0]
            infile = os.path.join(indir, f)
            gvcf = os.path.join(outdir, base.rsplit('_')[0])+'.g.vcf'

            '''
            --emit-ref-confidence/-ERC: ReferenceConfidenceModeMode for emitting
             reference confidence scores  Default value: NONE. Possible values:
             {NONE, BP_RESOLUTION, GVCF}

            ERC stands for Emit Reference Confidence.
            This means we want the program to estimate the probability of a
            given genotype being the reference base at each site in the file.
            GVCF mode tells the program to output every site in the genome
            (whether or not there appears to be a mutation).This is necessary
            for merging genotype information across multiple samples later.

            --sample-ploidy/-ploidy: Ploidy (number of chromosomes) per sample.
            For pooled data, set to (Number of samples in each pool * Sample Ploidy).
            Default value: 2.

            '''
            tool = 'HaplotypeCaller'
            java_options = repr('-Xmx4g -XX:ParallelGCThreads=1')
            erc = 'GVCF'

            if os.path.exists(gvcf):
                continue
            else:
                cmd = "java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T \
                {} -R {} -I {} --emitRefConfidence {} \
                -o {}".format(tool, reference, infile, erc, gvcf)
                try:
                    logging.info("calling variants in {}".format(os.path.basename(infile)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when calling variants in: {}".format(os.path.basename(infile)))
                    logging.error("Error that occurred: {}".format(e))
    return outdir

def genotype_gvcfs(indir, reference, outdir):
    """

    Gather multiple .gvcf files, by generating a single SNP file containing the
    genotypes at all sites with a sequence variant (mutation) in one or more
    individuals in our sample.

    Merges one or more HaplotypeCaller GVCF files into a single GVCF with
    appropriate annotations

    parameters:
        - indir (str): directory having every sample's gvcf file
        - reference (str): reference genome sequence in FASTA format
        - outdir (str): directory to store the output (a single/combined vcf file)

    returns:
        - outdir (str): directory to store the output (a single/combined vcf file)

    """

    tool = 'GenotypeGVCFs'
    java_options = repr('-Xmx4g -XX:ParallelGCThreads=1')

    gvcf_list = list()
    gvcf_dict = dict()
    for f in sorted(os.listdir(indir)):
        if f.endswith('.g.vcf'):
            f = os.path.join(indir, f)
            gvcf_list.append(f)

    for i, f in enumerate(gvcf_list):
        arg = '--variant {}'.format(f)
        if not i in gvcf_dict:
            gvcf_dict[i] = arg
    gvcfs = ' '.join(gvcf_dict.values())
    vcf_file = os.path.join(outdir, '{}samples.HC.vcf'.format(len(gvcf_list)))

    if os.path.exists(vcf_file):
        print("{} GVCF files already jointly genotyped into {}\
        ".format(str(len(gvcf_list)), os.path.basename(vcf_file)))
    else:
        cmd = ["java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T \
        {} -R {} {} \
        -o {}".format(tool, reference, gvcfs, vcf_file)]
        try:
            logging.info("performing joint genotyping on {} files\
            ".format(str(len(gvcf_list))))
            p = subprocess.Popen(cmd, shell=True)
        except Exception as e:
            logging.error("Error occurred when performing joint genotyping")
            logging.error("Error that occurred: {}".format(e))
    return outdir



############################## coverage ########################################
def assess_coverage(indir, reference, outdir):
    """
    Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library

    parameters:
        - indir (str): directory having BAM files
        - reference (str): reference genome sequence in FASTA format
        - outdir (str): directory to store DepthOfCoverage

    returns:
        - outdir (str): directory to store DepthOfCoverage
    """

    tool = 'DepthOfCoverage'
    partition = 'readgroup'

    for f in sorted(os.listdir(indir)):
        if f.endswith('fixmate.bam'):
            base = f.rsplit('.')[0]
            infile = os.path.join(indir, f)
            outfile = os.path.join(outdir, base.rsplit('_')[0])

            if os.path.exists(outfile):
                continue
            else:
                cmd = "java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T \
                {} -R {} -I {} -pt {} \
                -o {}".format(tool, reference, infile, partition, outfile)
                try:
                    logging.info("Assessing coverage in {}".format(os.path.basename(infile)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when assessing coverage in {}".format(os.path.basename(infile)))
                    logging.error("Error that occurred: {}".format(e))
    return outdir




##############################   data filtering ################################

def split_variants(indir, reference, outdir):
    """
    Create vcf files for only SNPs and only INDELS using SelectVariants
    (hard filtering)

    parameters:
        - indir(str): directory having combined vcf file from all samples
        - reference (str): reference genome sequence in FASTA format
        - outdir (str): directory to store the filtered variants
    returns:
        - outdir (str): directory to store the filtered variants
    """

    tool = 'SelectVariants'

    for f in sorted(os.listdir(indir)):
        if f.endswith('.vcf'):
            base = f.rsplit('.')[0]
            f = os.path.join(indir, f)
            snps_vcf = os.path.join(outdir, base)+'_raw_snps.vcf'
            indels_vcf = os.path.join(outdir, base)+'_raw_indels.vcf'

            if os.path.exists(snps_vcf) and os.path.exists(indels_vcf):
                print("SNP and INDEL variants already filtered, check {} & {}\
                ".format(os.path.basename(snps_vcf), os.path.basename(indels_vcf)))



            if not os.path.exists(snps_vcf):
                variant_type = 'SNP'
                try:
                    cmd = "java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T {} \
                    -R {} \
                    -V {} \
                    -selectType {} \
                    -o {}".format(tool, reference, f, variant_type, snps_vcf)
                    logging.info("hard-filtering {} in {}\
                    ".format(variant_type, os.path.basename(f)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when hard filtering SNP variants")
                    logging.error("Error that occurred: {}".format(e))

            if not os.path.exists(indels_vcf):
                variant_type = 'INDEL'
                try:
                    cmd = "java -jar /export/apps/gatk/3.7.0/GenomeAnalysisTK.jar -T {} \
                    -R {} \
                    -V {} \
                    -selectType {} \
                    -o {}".format(tool, reference, f, variant_type, indels_vcf)
                    logging.info("hard-filtering {} in {}\
                    ".format(variant_type, os.path.basename(f)))
                    p = subprocess.check_call(cmd, shell=True)
                except Exception as e:
                    logging.error("Error occurred when hard filtering SNP variants")
                    logging.error("Error that occurred: {}".format(e))
    return outdir

def main():
    # command line options
    parser=argparse.ArgumentParser()
    helpstr="""python snp_calling.py [options]"""
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('-d', '--datadir', required = True,
    help="directory having raw sequencing data (FASTQ format)", metavar='DIR')
    required_group.add_argument('-r', '--reference', required = True, help="reference/genome assembly file (FASTA format)")
    parser.add_argument('-t', '--threads', default=4, type=int, help="number of threads/cores")
    args=parser.parse_args()

    # check for raw data directory
    if args.datadir != None and os.path.exists(args.datadir) and os.path.isdir(args.datadir):
        DATA_DIR = args.datadir
        logging.info("Entering into the data directory {}".format(DATA_DIR))
    elif args.datadir != None and not os.path.exists(args.datadir):
        logging.error("Directory {0} doesn't exist".format(args.datadir))
        sys.exit(1)
    elif args.datadir != None and os.path.exists(args.datadir) and not os.path.isdir(args.datadir):
        logging.error("{} is not a directory".format(args.datadir))
    else:
        logging.error("please specify the directory having raw sequence data files\n")
        sys.exit(2)

    # check reference
    if args.reference != None and type(args.reference) == str:
        ref = args.reference
    else:
        logging.error("specify the reference sequence\n")

    # check number of threads
    if args.threads != None and type(args.threads) == int:
        thr = args.threads
    else:
        logging.error("specify the number of threads\n")


    # create the working directory
    user = getpass.getuser()
    project = os.path.basename(os.path.abspath(DATA_DIR))+'_output'
    WORKDIR = os.path.join(os.path.join('/var/scratch/', user), project)
    #WORKDIR = os.path.join(BASE_DIR, project)
    if not os.path.exists(WORKDIR):
        logging.info("Creating the working directory, all output will be stored here {}".format(WORKDIR))
        os.makedirs(WORKDIR)

    # create FASTQC directory
    FASTQC_DIR = os.path.join(WORKDIR, 'fastqc_output')
    if not os.path.exists(FASTQC_DIR):
        os.makedirs(FASTQC_DIR)
        logging.info("creating {} output directory".format(FASTQC_DIR))
    fastqc = run_fastqc(DATA_DIR, FASTQC_DIR, thr)

    # build reference index
    idx = build_genome_index(ref)

    # create BWA directory
    BWA_DIR = os.path.join(WORKDIR, 'bwa_output')
    if not os.path.exists(BWA_DIR):
        os.makedirs(BWA_DIR)
        logging.info("creating {} output directory".format(BWA_DIR))
    align = align_reads(fastqc, idx, BWA_DIR, thr)

    # sort and index SAM files
    sortindex = sortIndexsam(BWA_DIR)
    # collect alignments metrics
    metrics = collect_alignmentMetrics(BWA_DIR, ref)

    # create picard directory to store the deduplicated alignments
    PICARD_DIR = os.path.join(WORKDIR, 'picard_output')
    if not os.path.exists(PICARD_DIR):
        os.makedirs(PICARD_DIR)
        logging.info("creating {} output directory".format(PICARD_DIR))
    dedup = deduplicate_bam(BWA_DIR, PICARD_DIR, thr)

    # sort and index deduplicated BAM files
    sort_dedup = sortIndexsam(dedup)

    # validate the deduplicated and sorted BAM files
    validation = validate_bamfile(sort_dedup)

    # create reference index for GATK
    gatk_ref_idx = create_reference_index(ref)

    # define target intervals for realignment
    target_intervals = target_indel_intervals(PICARD_DIR, ref, thr)

    # realign
    realign = realign_indels(PICARD_DIR, ref, thr)

    # fix mate pair information
    fixmate = fix_matepair_info(PICARD_DIR)

    # SNP calling UnifiedGenotyper
    UG_DIR = os.path.join(WORKDIR, 'unified_genotyper_output')
    if not os.path.exists(UG_DIR):
        os.makedirs(UG_DIR)
        logging.info("creating {} output directory".format(UG_DIR))
    variants_ug = call_variants_UnifiedGenotyper(PICARD_DIR, ref, UG_DIR, thr)

    # single SNP calling with HaplotypeCaller
    # HC_DIR = os.path.join(WORKDIR, 'haplotypecaller_output')
    # if not os.path.exists(HC_DIR):
    #     os.makedirs(HC_DIR)
    #     logging.info("creating {} output directory".format(HC_DIR))
    # variants_hc = call_variants_HaplotypeCaller(PICARD_DIR, ref, HC_DIR)
    #
    #
    # # create GATK directory and gather variants
    # VCF_DIR = os.path.join(WORKDIR, 'vcf_output')
    # if not os.path.exists(VCF_DIR):
    #     os.makedirs(VCF_DIR)
    #     logging.info("creating {} output directory".format(VCF_DIR))
    # vcf = genotype_gvcfs(HC_DIR, ref, VCF_DIR)
    #
    #
    # create coverage output and assess depth of coverage
    DP_DIR = os.path.join(WORKDIR, 'coverage_output')
    if not os.path.exists(DP_DIR):
        os.makedirs(DP_DIR)
        logging.info("creating {} output directory".format(DP_DIR))
    depth = assess_coverage(PICARD_DIR, ref, DP_DIR)


    # create directory and filter variants
    FILTERED_DIR = os.path.join(WORKDIR, 'filtered_variants_output')
    if not os.path.exists(FILTERED_DIR):
        os.makedirs(FILTERED_DIR)
        logging.info("creating {} output directory".format(FILTERED_DIR))
    filter_variants = split_variants(variants_ug, ref , FILTERED_DIR)




    # set end time format
    end_time = datetime.now()
    e = end_time.strftime(format)
    tdelta = end_time - start_time



    # format the time delta object to human readable form
    d = dict(days=tdelta.days)
    d['hrs'], rem = divmod(tdelta.seconds, 3600)
    d['min'], d['sec'] = divmod(rem, 60)

    if d['min'] is 0:
        fmt = '{sec} sec'
    elif d['hrs'] is 0:
        fmt = '{min} min {sec} sec'
    elif d['days'] is 0:
        fmt = '{hrs} hr(s) {min} min {sec} sec'
    else:
        fmt = '{days} day(s) {hrs} hr(s) {min} min {sec} sec'


    logging.info("\nCompleted analysis on: "+e+ "\n" )
    logging.info("[ALL done] Runtime: " +'\t'+fmt.format(**d))
if __name__ == '__main__':
    main()
