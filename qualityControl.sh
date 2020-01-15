#!/usr/bin/bash

#SBATCH -p batch
#SBATCH -w compute04
#SBATCH -J QC
#SBATCH -n 4
#SBATCH -o slurm_out
#SBATCH -e slurm_errors

#loading fastqc
module load fastqc/0.11.7

BASEDIR=/var/scratch/jb/Geoffrey_rice_blast/Data
OUTDIR=/var/scratch/mwanjiku

mkdir -p ${OUTDIR}/FastQC_reports

for i in 1- 11 12 12- 13 14 15 15- 16 16- 17 18 21- 22 22- 23- 24 25 27 29 32 32- 34 36 38 39- 4- 41 43 47 49 5 5- 52 54 55 56 58 59 6 60 61 64 65 69 7 7- 70 71 9
do
   echo -e "\n\nIsolate ${i}\n"
   cd ${BASEDIR}/Isolate${i}
   fastqc ${BASEDIR}/Isolate${i}/*.fq.gz -o ${OUTDIR}/FastQC_reports

done


