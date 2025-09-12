#! /bin/bash
#SBATCH -p barc
#SBATCH -n 12
#SBATCH --mem=100G
#SBATCH -t 2:00:00

module load python/3.12
source .venv/bin/activate

python model/genotype_loader.py