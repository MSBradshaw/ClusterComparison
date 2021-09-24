#!/bin/bash

#SBATCH --job-name snowball_booster
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --partition shas
#SBATCH --qos normal
#SBATCH --time=24:00:00
#SBATCH --chdir=/scratch/summit/cgibbs10@colostate.edu/projs/ClusterComparison/
#SBATCH -o msgs.output
#SBATCH -e msgs.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=connor.gibbs@colostate.edu

module purge
module load anaconda

Rscript AnalyzeResults/snowball_booster.R