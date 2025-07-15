
def chrs = (1..22).collect { it.toString() } + ['X']
//def chrs = (22).collect { it.toString() } 

def files = chrs.collect { chr -> 
    def bgen_file = file("${params.dir}/ukb_imp_chr${chr}_bct_vars.bgen") 
    def sample_file = file("${params.dir}/ukb_imp_chr${chr}_bct_vars.sample") 
    return[bgen_file, sample_file]
}

def grouped_files = Channel.from(files)

process bgenToVCF {

    input:
    tuple path(bgen_file), path(sample_file)

    output:
    file "${bgen_file.baseName}.vcf"

    script:
    """
    module load plink    
    plink2 --bgen ${bgen_file} ref-first --sample ${sample_file} --export vcf --memory ${params.max_memory} --threads ${params.threads} --out ${bgen_file.baseName}
    """

}

process cleanVCF {
    input: 
    path vcf_file

    output:
    file "${vcf_file.baseName}_clean.vcf"
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    awk 'BEGIN { OFS="\t" }
    {
        if(/^##/ || (NR==1 && \$0 ~ /^#/)) {
        print \$0;
        } else {
        print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8;
        }
    }' "${vcf_file}" > "${vcf_file.baseName}_clean.vcf"
    """
}

process snpEff {
    input:
    path vcf_file

    output:
    file "${vcf_file.baseName}" //summary table
    file "${vcf_file.baseName}.genes.txt"
    file "${vcf_file.baseName}.ann.vcf", emit: annotated_vcf
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load snpeff
    snpEff -v ${params.reference} ${vcf_file} -dataDir ${projectDir}/../${params.data_directory} -csvStats ${vcf_file.baseName} > ${vcf_file.baseName}.ann.vcf
    """
}

process compileVariants {

    input:
    path vcf_files

    output:
    file "bct_variant_annotated.csv"
    publishDir "${params.outdir}/data", mode: 'copy'

    """
    #!${projectDir}/.venv/bin/python3
    
    import sys
    import os
    import glob
    sys.path.append('${projectDir}')
    
    import src.vcf as vcf
    import pandas as pd

    # Process all annotated VCF files
    dfs = []
    vcf_files = glob.glob('*.ann.vcf')
    
    for vcf_file in vcf_files:
        print(f'Processing {vcf_file}...')
        variant_dict = vcf.readVCF(vcf_file)
        dfs.append(pd.DataFrame.from_dict(variant_dict))

    # Concatenate all dataframes
    if dfs:
        df = pd.concat(dfs, ignore_index=True)
        df.to_csv('bct_variant_annotated.csv', index=False)
        
        unique_variants = len(df['id'].unique())
        print(f'{unique_variants} unique variants after compiling')
    else:
        print('No annotated VCF files found')
    """
}


workflow {

    if (params.skip_vcf) {
        def vcf_files = Channel.from(params.vcf_files) 
        def annotated_vcfs = snpEff(vcf_files)
    } else {
        def annotated_vcfs = bgenToVCF(grouped_files) |
        cleanVCF |
        snpEff
    }

    annotated_vcfs.annotated_vcf.collect() |
    compileVariants
}