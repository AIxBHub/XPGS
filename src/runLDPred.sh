#!/bin/bash
#SBATCH -p barc
#SBATCH --constraint=rhel9
#SBATCH --mem=200G
#SBATCH -n 12
#SBATCH -t 3:00:00

module load apptainer

apptainer exec container/bioc_glm_bioc320_R442_amd64.sif Rscript src/ldpred2.R