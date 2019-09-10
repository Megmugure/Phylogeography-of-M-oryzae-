#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -w compute04
#SBATCH -J spades_jb
#SBATCH -n 4
#SBATCH -o /var/scratch/jb/Geoffrey_rice_blast/Logs/spades_job_log
#SBATCH -e /var/scratch/jb/Geoffrey_rice_blast/Logs/spades_job_errlog

# load the spades module
module load spades/3.11.1

BASEDIR=/var/scratch/jb/Geoffrey_rice_blast/Data # this is on compute04
REL_OUTDIR=Spades_output # relative to the Isolatedir


for ISOLATE in 1- 11 12 12- 13 14 15 15- 16 16- 17 18 21- 22 22- 23- 24 25 27 29 32 32- 34 36 38 39- 4- 41 43 47 49 5 5- 52 54 55 56 58 59 6 60 61 64 65 69 7- 70 71 9 # all except 7, that is already done
do
  ISOLATEDIR=${BASEDIR}/Isolate${ISOLATE}/
  cd ${ISOLATEDIR}

  # making directory for output
  mkdir -p ${REL_OUTDIR} # no error if existing

  # spades to assemble all the reads
  spades.py --careful --phred-offset 33 --threads 4 \
          --pe1-1 *_1.fq.gz \
          --pe1-2 *_2.fq.gz \
          -o ${REL_OUTDIR}
done
