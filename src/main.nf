
def chrs = (1..22).collect { it.toString() } + ['X']
//def chrs = (22).collect { it.toString() } 

def files = chrs.collect { chr -> 
    def bgen_file = file("${params.dir}/ukb_imp_chr${chr}_bct_vars.bgen") 
    def sample_file = file("${params.dir}/ukb_imp_chr${chr}_bct_vars.sample") 
    return[bgen_file, sample_file]
}

def grouped_files = Channel.from(files)

process extractPGSVariants {
    
    input:
    val pgs_ids
    
    output:
    path "pgs_variants.txt"
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    # PGS IDs to process
    pgs_ids = "${pgs_ids}".split(',')
    all_variants = set()
    
    # Base URL for PGS scoring files
    base_url = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{}/ScoringFiles/Harmonized/{}_hmPOS_GRCh38.txt.gz"
    
    dtype_dict = {
        'rsID': str,
        'chr_name': str, 
        'chr_position': 'Int64',
        'hm_chr': str,
        'hm_pos': 'Int64',
        'hm_rsID': str,
    }
    
    for pgs_id in pgs_ids:
        pgs_id = pgs_id.strip()
        url = base_url.format(pgs_id, pgs_id)
        print(f"Downloading variants from {pgs_id}...")
        
        try:
            df = pd.read_csv(url, 
                           compression='gzip',
                           sep='\t', 
                           comment='#',
                           low_memory=False,
                           dtype=dtype_dict,
                           na_values=['NA', 'na', ''])
            
            # Extract rsIDs, removing any NaN values
            rsids = df['rsID'].dropna().unique()
            all_variants.update(rsids)
            print(f"Added {len(rsids)} variants from {pgs_id}")
            
        except Exception as e:
            print(f"Warning: Could not download {pgs_id}: {e}")
    
    # Write variant list for plink2
    print(f"Total unique variants: {len(all_variants)}")
    with open('pgs_variants.txt', 'w') as f:
        for variant in sorted(all_variants):
            f.write(f"{variant}\\n")
    """
}

process bgenToVCF {

    input:
    tuple path(bgen_file), path(sample_file), path(variant_list)

    output:
    path("${bgen_file.baseName}.vcf"), optional: true

    script:
    """
    module load plink    
    plink2 --bgen ${bgen_file} ref-first --sample ${sample_file} \
      --extract ${variant_list} \
      --export vcf \
      --memory ${params.max_memory} --threads ${task.cpus} \
      --out ${bgen_file.baseName}
    """

}

process sortVCF {
    input:
    path vcf_file

    output:
    file "${vcf_file.baseName}_sorted.vcf.gz"
//    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load bcftools
    bcftools sort ${vcf_file} -Oz -o ${vcf_file.baseName}_sorted.vcf.gz
    """
}

//process indexVCF {
//    input:
//    path vcf_file
//
//    output:
//    file "${vcf_file}.tbi"
//
//    script:
//    """
//    module load bcftools
//    bcftools index --tbi ${vcf_file}
//    """
//}
////temp -- should incorporate compression with sort step
//process compressVCF {
//    stageInMode 'link' //pigz won't work with symbolic links
//    
//    input:
//    path vcf_file
//
//    output:
//    file "${vcf_file}.gz"
//
//    script:
//    """
//    pigz -p ${task.cpus} ${vcf_file}
//    """
//}

process mergeVCF {
    input:
    path vcf_files
//    path indexed_vcfs

    output:
    file "ukb_imp_bct_vars_merged.vcf"
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load bcftools
    bcftools concat ${vcf_files} -Ov -o ukb_imp_bct_vars_merged.vcf
    """
}

process vcfToBgen {
    input:
    path vcf_file

    output:
    file "${vcf_file.simpleName}.bed"
    file "${vcf_file.simpleName}.bim"
    file "${vcf_file.simpleName}.fam"
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load plink/2.00
    plink2 --vcf ${vcf_file} --make-bed --split-par hg38 --threads ${params.threads} --out ${vcf_file.simpleName}
    """
}

process cleanVCF {
    input: 
    path vcf_file

    output:
    file "${vcf_file.simpleName}_clean.vcf"
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
    }' "${vcf_file}" > "${vcf_file.simpleName}_clean.vcf"
    """
}

process snpEff {
    input:
    path vcf_file

    output:
    file "${vcf_file.simpleName}" //summary table
    file "${vcf_file.simpleName}.genes.txt"
    file "${vcf_file.simpleName}.ann.vcf"
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load snpeff
    snpEff -v ${params.reference} ${vcf_file} -dataDir ${projectDir}/../${params.data_directory} -csvStats ${vcf_file.simpleName} > ${vcf_file.simpleName}.ann.vcf
    """
}

workflow {

//    if (params.skip_vcf) {
//        def vcf_files = Channel.from(params.vcf_files) 
//        snpEff(vcf_files)
//    } else {
        
    def variant_list = extractPGSVariants(params.pgs_ids)
    // combine input bgen with the exctracted variant file
    // then map the plink files and variants to create a tuple input for bgenToVCF
    def merged_vcfs = grouped_files.combine(variant_list) |
        map { bgen, sample, variants -> tuple(bgen, sample, variants) } |
        bgenToVCF | 
        sortVCF |
        collect |
        mergeVCF
    
    //def indexed_vcfs = indexVCF(sorted_vcfs) 
    
    //def merged_vcfs = mergeVCF(sorted_vcfs.collect(), indexed_vcfs.collect())
    // steps to make annotation files
    cleanVCF(merged_vcfs) | snpEff
    
    // converty back out to binary 
    vcfToBgen(merged_vcfs)
}