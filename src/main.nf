
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

process mergeVCF {
    input:
    path vcf_files

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
    plink2 --vcf ${vcf_file} --make-bed --split-par hg38 --memory 500000 --threads ${task.cpus} --out ${vcf_file.simpleName}
    """
}

// TODO -- fix this; there is a param in bcftools to drop genotypes with subset -G (after subsetting if -s is set)
process clean_vcf {
    input: 
    path vcf_file

    output:
    file "${vcf_file.simpleName}_clean.vcf"
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load bcftools
    bcftools view -G ${vcf_file} -o ${vcf_file.simpleName}_clean.vcf 
    """
}

process snpEff {
    input:
    path vcf_file

    output:
    tuple path("${vcf_file.simpleName}"), path("${vcf_file.simpleName}.genes.txt"), path("${vcf_file.simpleName}.ann.vcf")
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load snpeff
    snpEff -v ${params.reference} ${vcf_file} -dataDir ${projectDir}/../${params.data_directory} -csvStats ${vcf_file.simpleName} > ${vcf_file.simpleName}.ann.vcf
    """
}

process generateAnnotationCSV {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path annotated_vcf

    output:
    path "${annotated_vcf.simpleName}_annotations.csv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import vcfpy
    import os

    def readVCF(vcf):
        dict = {
            'chrom' : [],
            'pos' : [],
            'id' : [],
            'ann' : [],
            'gene' : [],
            'biotype' : []
        }

        vcf_reader = vcfpy.Reader.from_path(vcf)

        for line in vcf_reader:
            ## INFO -- dictionary with ANN as key to list of '|' delimited strings, each annotation is an element in the list
            annotations = line.INFO.get('ANN',[]) 

            ## extract variants that have annotations to protein coding genes
            annotations = [a for a in annotations if 'protein_coding' in a or 'intergenic_region' in a] 

            for anno in annotations: ## then loop the annotations
     
                fields = anno.split('|')

                dict['chrom'].append(line.CHROM)
                dict['pos'].append(line.POS)
                dict['id'].append(line.ID[0])
                
                dict['ann'].append(fields[1])
                dict['gene'].append(fields[4])
                dict['biotype'].append(fields[7])

        return(dict)
    
    # Read VCF annotations using the readVCF function
    variant_dict = readVCF('${annotated_vcf}')
    
    # Convert to DataFrame and save as CSV
    df = pd.DataFrame.from_dict(variant_dict).drop_duplicates()

    effects = "${params.effect}".split(',')
    print(effects)
    # grab variants of interest based on specific annotation(s) from params
    df_VOI = df[df['ann'].isin(effects)]
    df_VOI.to_csv('${annotated_vcf.simpleName}_effects.csv', index = False)
    # save out all the annotations
    df.to_csv('${annotated_vcf.simpleName}_annotations.csv', index=False)
    
    """
}

process subset_vcf {
    input:
    path annotated_vcf
    val effect_expression

    output:
    path ("${annotated_vcf.simpleName}_effects.vcf")
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load snpeff
    echo "${effect_expression}" > effect_expression.txt
    snpSift filter -e effect_expression.txt ${annotated_vcf} > ${annotated_vcf.simpleName}_effects.vcf
    """
}

// Function to build effect expression from a list of effects
def buildEffectExpression(effects) {
    return effects.collect { effect -> "ANN[*].EFFECT has '${effect}'" }.join(' | ')
}

workflow {
    // Create a value channel for the effect expression
    def effectsList = params.effect instanceof List ? params.effect : params.effect.split(',').collect { it.trim() }
    def effectExpression = buildEffectExpression(effectsList)
    def effectExpressionCh = Channel.value(effectExpression)

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
    def snpeff_results = snpEff(merged_vcfs)
    
    clean_vcf(snpeff_results.map { stats, genes, vcf -> vcf }) | // drop genotypes from the annotated vcf
        generateAnnotationCSV

    // Use the effect expression channel with the subset_vcf process
    subset_vcf(snpeff_results.map { stats, genes, vcf -> vcf }, effectExpressionCh) // subset to specified effect, retaining genotypes
    // converty back out to binary 
    vcfToBgen(merged_vcfs)
}