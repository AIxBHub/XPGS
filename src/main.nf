
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
 //   publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load plink    
    plink2 --bgen ${bgen_file} ref-first --sample ${sample_file} --export vcf --memory ${params.max_memory} --threads ${params.threads} --out ${bgen_file.baseName}
    """

}

process sortVCF {
    input:
    path vcf_file

    output:
    file "${vcf_file.baseName}_sorted.vcf"
//    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load bcftools
    bcftools sort ${vcf_file} -o ${vcf_file.baseName}_sorted.vcf
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
    bcftools concat ${vcf_files} -o ukb_imp_bct_vars_merged.vcf
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
    plink2 --vcf ${vcf_file} --make-bed --split-par --threads ${params.threads} --out ${vcf_file.simpleName}
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

    if (params.skip_vcf) {
        def vcf_files = Channel.from(params.vcf_files) 
        snpEff(vcf_files)
    } else {
        
        def merged_vcf_files = bgenToVCF(grouped_files) |
        sortVCF |
        collect |
        mergeVCF

        cleanVCF(merged_vcf_files) |
        snpEff

        vcfToBgen(merged_vcf_files)
    }
}