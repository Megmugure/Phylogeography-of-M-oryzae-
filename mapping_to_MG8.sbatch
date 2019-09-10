#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -w compute04
#SBATCH -J mapping
#SBATCH -n 4

# load the necessary modules
module load bowtie2/2.3.4.1
module load samtools/1.8

BASEDIR=/var/scratch/jb/Geoffrey_rice_blast/Data # this is on compute04
MG8REFDIR=/var/scratch/jb/Geoffrey_rice_blast/Magnaporthe_references/MG8_GCA_000002495.2
MG8REF_FASTA=${MG8REFDIR}/Magnaporthe_oryzae.MG8.dna.toplevel.fa

# first, make the index for Bowtie2 (once-off)
# cd ${MG8REFDIR}
# bowtie2-build --threads 2 ${MG8REF_FASTA} MG8_index

for i in 1- 11 12 12- 13 14 15 15- 16 16- 17 18 21- 22 22- 23- 24 25 27 29 32 32- 34 36 38 39- 4- 41 43 47 49 5 5- 52 54 55 56 58 59 6 60 61 64 65 69 7 7- 70 71 9
do
  if [ -e "${BASEDIR}/Isolate${i}/bowtie.sam" ]
  then
    continue
  fi

  echo -e "\n\nIsolate ${i}...\n"
  cd ${BASEDIR}/Isolate${i}

  bowtie2 --phred33 --threads 4 -x ${MG8REFDIR}/MG8_index -1 *_1.fq.gz -2 *_2.fq.gz -S bowtie.sam 

  # the following is some code to use with unclean reads, to trim them first.
  # bowtie2 to align PAIRED survivor reads from trimmomatic
  # bowtie2 --threads 4 --un-conc sample${i}_P_unmapped.fastq --no-unal --no-discordant -x $PIGREF -1 trimmed_sample${i}_1P -2 trimmed_sample${i}_2P -S sample${i}_P_mapped.sam &>> bowtie2.log
  # option --no-discordant to prevent search of discordant alignments altogether, saving time.
  # bowtie2 to align UNPAIRED survivor reads from trimmomatic
  # bowtie2 --threads 4 --un sample${i}_U_unmapped.fastq --no-unal -x $PIGREF -U trimmed_sample${i}_1U,trimmed_sample${i}_2U -S sample${i}_U_mapped.sam &>> bowtie2.log

  # SAMTOOLS to create BAM file
  samtools view -@ 3 -bT ${MG8REF_FASTA} -o bowtie_isolate${i}.bam bowtie.sam # -@ 3 means 4 threads in total

  samtools sort -@ 4 bowtie_isolate${i}.bam > bowtie_isolate${i}_sorted.bam # here for samtools sort, -@ is followed by the total # of cores to use

  samtools index bowtie_isolate${i}_sorted.bam # the index file is necessary to tools e.g. Tablet
  # samtools index doesn't seem to be multi-threaded

done
