#! /usr/bin/env bash
#SBATCH -p batch
#SBATCH -w compute05
#SBATCH -J converter2
#SBATCH -n 8

#SBATCH -o /home/mwanjiku/R_analysis/slurm_out
#SBATCH -e /home/mwanjiku/R_analysis/slurm_errors

# loading the R module
module load R/3.6

# create an output directory
#mkdir -p /home/mwanjiku/R_analysis/R_output
cd /home/mwanjiku/R_analysis/R_output


# running the R script
Rscript /home/mwanjiku/scripts/adegenet2.R

