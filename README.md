## Blood Cell Trait Variant Processing

Workflow to...
1. Find pgs models from pgs_catalog for all blood cell traits 
2. Filter models for performance and European ancestry
3. Generate a union of variant IDs to get imputed genotypes

## Nextflow Pipeline for variant subsetting and annotation

Pipeline takes the imputed genotypes (per chromosome plink files)...
1. Pull pgs score file with pgs model id and extract target variants in the process
2. Convert bgen to vcf with `plink2` and extracts 
3. sorts and compress vcf files (per chromosome)
4. merges per chrosome vcf (`bcftools concat`)
5. Annotate merged vcf with `snpEff`
6. Converts merged vcf file back to binary -- for ldpred
7. Generates csv lookup file with processed annotations

# Running

## Run PGS processing for union variants

Union variants for blood cell traits (bct) were identified by running `pgs_utils.py`.
Trait ids are specified in a csv files `sheets/traits.csv`. 
Specify sheet and outdirectory in `__main__` of `pgs_utils.py`

```
module load python/3.12
source .venv/bin/activate

python src/pgs_utils.py
```

## Run Nextflow 

Imputed genotypes were processed with nextflow. 

`config/nfparams.yaml` specifies path to bgen files and reference for `snpEff`. 
- specify model id(s) with `pgs_ids`
- specify outdirectory with `outdir`
- specify path to binary plink files `dir`
- specify reference genome for `snpEff` with `reference`

nextflow configurations are in `config/nextflow.config` for specific process configs. 

Apptainer cache is set default to `.cach` in the working directory.


Then run with...
```
sbatch src/runNF.sh
```
# snpEff reference info

config file: `/nas/{CLUSTERNAME}/rhel8/apps/snpeff/5.2c/snpEff/snpEff.config`

# Notes

#TODO add instructions for `nn_id_from_snpeff.py`

summary stats
https://pmc.ncbi.nlm.nih.gov/articles/PMC7482360/
https://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/