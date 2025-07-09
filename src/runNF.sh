#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --mem=8G
#SBATCH -p barc
#SBATCH --constraint=rhel9

module load nextflow

nextflow run src/main.nf -resume -params-file config/nfparams.yaml -c config/nextflow.config