# Blood Cell Trait Variant Processing

Workflow to...
1. Find pgs models from pgs_catalog for all blood cell traits 
2. Filter models for performance and European ancestry
3. Generate a union of variant IDs to get imputed genotypes

With imputed genotypes...
1. Convert bgen to vcf with `plink2`
2. Annotate vcf with `snpEff`
3. Standardize gene ids in annotation with `NodeNorm`
4. Subset annotated variants for specific pgs models -- *should have been used to make the variants union!*

# Running

Union variants for blood cell traits (bct) were identified by running `pgs_utils.py`.
Trait ids are specified in a csv files `sheets/traits.csv`. 
Specify sheet and outdirectory in `__main__` of `pgs_utils.py`

```
module load python/3.12
source .venv/bin/activate

python src/pgs_utils.py
```

Imputed genotypes were processed with nextflow. 
`config/nfparams.yaml` specifies path to bgen files and reference for `snpEff`. 
There are a few other params there, including an option to skip vcf conversion in case you just want to reannotate with a different reference. 

```
sbatch src/runNF.sh
```

#TODO add instructions for `nn_id_from_snpeff.py`

Annotated variants can be subset to specific pgs models like so.
```
import src.vcf as vcf
subset = vcf.subsetVariants({PGS_MODEL_ID}, {PATH_TO_ANNOTATED_VARIANTS})
subset.to_csv({SUBSET_OUTPUT_FILE})
```
#TODO clean up `vcf.py` 

# snpEff reference info

config file: `/nas/{CLUSTERNAME}/rhel8/apps/snpeff/5.2c/snpEff/snpEff.config`
