#!/usr/bin/bash

#SBATCH -p batch
#SBATCH -w compute05
#SBATCH -J filtering
#SBATCH -n 4
#SBATCH -o /home/mwanjiku/filtering_out/filtering_out.log
#SBATCH -e /home/mwanjiku/filtering_out/filtering_errors.log


BASEDIR=/home/mwanjiku/variantCalling_out
OUTDIR=/home/mwanjiku/filtering_out

#load the necessary modules
module load vcftools/0.1.15

for i in 11 12 12- 13 14 15 15- 16- 17 1- 18 21- 22 22- 23- 24 25 27 29 32 34 36 38 39- 4- 41 43 47 49 5 5- 52 54 55 56 59 6 60 69 70 71 9


do

   vcftools --vcf ${BASEDIR}/bowtie_isolate${i}.vcf --remove-indels --recode --recode-INFO-all --out ${OUTDIR}/isolate${i}
 
   

done
