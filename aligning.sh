#!/usr/bin/bash

#SBATCH -p batch
#SBATCH -w compute04
#SBATCH -J mapping
#SBATCH -n 8
#SBATCH -o slurm_out
#SBATCH -e slurm_errors


#loading bowtie and samtools
module load bowtie2/2.3.4.1
module load samtools/1.8

BASEDIR=/home/mwanjiku/missing_reads/blast_reads
OUTDIR=/home/mwanjiku
MG8REFDIR=/home/mwanjiku/refGenome
MG8REF_FASTA=${MG8REFDIR}/Magnaporthe_oryzae.MG8.dna.toplevel.fa

building an index of the reference genome for bowtie2
cd ${MG8REFDIR}
bowtie2-build --threads 8 ${MG8REF_FASTA} refIndex

#mapping the reads to the indexed reference genome
#mkdir -p ${OUTDIR}/aligning_out
for i in 16 58 61 64 65 7 7-
do

  echo -e "\n\nIsolate ${i}\n"
  cd ${BASEDIR}/Isolate${i}

  /export/apps/bowtie2/2.3.4.1/bowtie2-align-s --wrapper basic-0 --phred33 --threads 8 -x ${MG8REFDIR}/refIndex -1 ${BASEDIR}/Isolate${i}/*_1.fq.gz -2 ${BASEDIR}/Isolate${i}/*_2.fq.gz -S ${OUTDIR}/aligning_out/bowtie${i}.sam  ##use of samtools to create a BAM file


  samtools view -@ 8 -bT ${MG8REF_FASTA} -o ${OUTDIR}/aligning_out/bowtie_isolate${i}.bam ${OUTDIR}/aligning_out/bowtie${i}.sam
  samtools sort -@ 8 -o ${OUTDIR}/aligning_out/bowtie_isolate${i}_sorted.bam ${OUTDIR}/aligning_out/bowtie_isolate${i}.bam
  samtools index ${OUTDIR}/aligning_out/bowtie_isolate${i}_sorted.bam
done

