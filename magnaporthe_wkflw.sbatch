#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -w compute04
#SBATCH -J magnaporthe
#SBATCH -n 4

# load the necessary modules
module load bowtie2/2.3.4.1
module load fastqc/0.11.7
module load trimmomatic/0.38
module load samtools/1.8

BASEDIR=/var/scratch/jb/Geoffrey_rice_blast # this is on compute04
MG8REFDIR=/var/scratch/jb/Geoffrey_rice_blast/Magnaporthe_references/MG8_GCA_000002495.2
MG8REF_FASTA=${MG8REFDIR}/Magnaporthe_oryzae.MG8.dna.toplevel.fa

# first, make the index for Bowtie2 (once-off)
# cd ${MG8REFDIR}
# bowtie2-build --threads 2 ${MG8REF_FASTA} MG8_index

# mkdir -p ${BASEDIR}/FastQC_reports

for i in 7 32- 58 61 64 65 # the isolates of interest
# for i in 7 7- 32- 58 61 64 65 # the isolates of interest
# for i in 7
do
  echo -e "\n\nIsolate ${i}...\n"
  # FASTQC
  # fastqc --threads 2 -o ${BASEDIR}/FastQC_reports -d ${BASEDIR} ${BASEDIR}/Data/Isolate${i}/*.fq.gz

  # TRIMMOMATIC: light trimming
  # trimmomatic PE -threads 4 -phred33 -basein xxxxx -baseout xxxxx ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 # these trimmers are the default ones from http://www.usadellab.org/cms/?page=trimmomatic
  # in our case here we decide not to trim anything, since BGI has taken care of adapter trimming and low-quality reads removal

  # BOWTIE2
  cd ${BASEDIR}/Data/Isolate${i}
  # bowtie2 --phred33 --threads 4 -x ${MG8REFDIR}/MG8_index -1 *_1.fq.gz -2 *_2.fq.gz -S bowtie.sam 

  # the following is some code to use with unclean reads, to trim them first.
  # bowtie2 to align PAIRED survivor reads from trimmomatic
  # bowtie2 --threads 4 --un-conc sample${i}_P_unmapped.fastq --no-unal --no-discordant -x $PIGREF -1 trimmed_sample${i}_1P -2 trimmed_sample${i}_2P -S sample${i}_P_mapped.sam &>> bowtie2.log
  # option --no-discordant to prevent search of discordant alignments altogether, saving time.
  # bowtie2 to align UNPAIRED survivor reads from trimmomatic
  # bowtie2 --threads 4 --un sample${i}_U_unmapped.fastq --no-unal -x $PIGREF -U trimmed_sample${i}_1U,trimmed_sample${i}_2U -S sample${i}_U_mapped.sam &>> bowtie2.log

  # SAMTOOLS to create BAM file
  # samtools view -bT ${MG8REF_FASTA} -o bowtie_isolate${i}.bam bowtie.sam
  samtools sort bowtie_isolate${i}.bam > bowtie_isolate${i}_sorted.bam
  samtools index bowtie_isolate${i}_sorted.bam # the index file is necessary to tools e.g. Tablet

done
