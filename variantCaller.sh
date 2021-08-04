#!/usr/bin/bash

#SBATCH -p batch
#SBATCH -w compute04
#SBATCH -J bcf
#SBATCH -n 8
#SBATCH -o variant_out.log
#SBATCH -e variant_errors.log
THREADS=8


(( BCF_THREADS = THREADS-1 ))

BASEDIR=/home/mwanjiku/aligning_out
MG8REFDIR=/home/mwanjiku/refGenome/Magnaporthe_oryzae.MG8.dna.toplevel.fa
OUTDIR=/home/mwanjiku/variantCalling_out

#loading bcftools and igvtools

module load bcftools/1.8
module load igvtools/2.3.98

#running the following commands on all the data

# running the commands on only the 42 available isolates
# for i in 12 12- 13 14 15 15- 16- 17 18 21- 22 22- 23- 24 25 27 29 32 34 36 38 39- 4- 41 43 47 49 5 5- 52 54 55 56 59 6 60 69 70 71 9
# running the commands on the remaining 7 isolates
for i in 16 58 61 64 65 7 7-

do
  echo -e "\n\nIsolate ${i}\n"

  bcftools mpileup --threads ${BCF_THREADS} -Ob -o ${OUTDIR}/bowtie_isolate${i}.bcf -f ${MG8REFDIR} ${BASEDIR}/bowtie_isolate${i}_sorted.bam

  bcftools call -vmO z -o ${OUTDIR}/bowtie_isolate${i}.vcf.gz ${OUTDIR}/bowtie_isolate${i}.bcf

  igvtools index ${OUTDIR}/bowtie_isolate${i}.vcf.gz
 

done
