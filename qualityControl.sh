#!/usr/bin/bash

#SBATCH -p batch
#SBATCH -w compute04
#SBATCH -J QC
#SBATCH -n 4
#SBATCH -o slurm_out
#SBATCH -e slurm_errors
THREADS=4


#loading fastqc
module load fastqc/0.11.7

((QC_THREADS=THREADS-1))

BASEDIR=/home/mwanjiku/missing_reads/blast_reads
OUTDIR=/home/mwanjiku

mkdir -p ${OUTDIR}/FastQC_reports

for i in 16 58 61 64 65 7 7-
do
   echo -e "\n\nIsolate ${i}\n"
   cd ${BASEDIR}/Isolate${i}
   fastqc ${BASEDIR}/Isolate${i}/*.fq.gz --threads ${QC_THREADS} -o ${OUTDIR}/FastQC_reports

done


