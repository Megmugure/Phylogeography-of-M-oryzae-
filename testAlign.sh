#!/usr/bin/bash

#SBATCH -p batch
#SBATCH -w compute04
#SBATCH -J mapping
#SBATCH -n 8
#SBATCH -o aligning_out.log
#SBATCH -e aligning_errors.log

#loading bowtie and samtools
module load bowtie2/2.3.4.1
module load samtools/1.8

BASEDIR=/var/scratch/jb/Geoffrey_rice_blast/Data
OUTDIR=/var/scratch/mwanjiku
MG8REFDIR=/home/mwanjiku/refGenome
MG8REF_FASTA=${MG8REFDIR}/Magnaporthe_oryzae.MG8.dna.toplevel.fa

#building an index of the reference genome for bowtie2
#cd ${MG8REFDIR}
#bowtie2-build --threads 4 ${MG8REF_FASTA} refIndex

#mapping the reads to the indexed reference genome
mkdir -p ${OUTDIR}/testAlign_out
for i in  39-
do
  if [ -e "${OUTDIR}/testAlign_out/bowtie${i}.sam" ]
  then
    continue
  fi

    echo -e "\n\nIsolate ${i}\n"
    cd ${BASEDIR}/Isolate${i}

  /export/apps/bowtie2/2.3.4.1/bowtie2-align-s --wrapper basic-0 --phred33 --threads 4 -x ${MG8REFDIR}/refIndex -1 ${BASEDIR}/Isolate${i}/*_1.fq.gz -2 ${BASEDIR}/Isolate${i}/*_2.fq.gz -S ${OUTDIR}/testAlign_out/bowtie${i}.sam  ##use of samtools to create a BAM file


  samtools view -@ 3 -bT ${MG8REF_FASTA} -o ${OUTDIR}/testAlign_out/bowtie_isolate${i}.bam ${OUTDIR}/testAlign_out/bowtie${i}.sam
  samtools sort -@ 4 ${OUTDIR}/testAlign_out/bowtie_isolate${i}.bam > ${OUTDIR}/testAlign_out/bowtie_isolate${i}_sorted.bam
  samtools index ${OUTDIR}/testAlign_out/bowtie_isolate${i}_sorted.bam ${OUTDIR}/testAlign_out/bowtie_isolate${i}_sorted.bam.bai
done


