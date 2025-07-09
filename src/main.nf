
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
    file "${vcf_file.baseName}.ann.vcf"
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    module load snpeff
    snpEff -v ${params.reference} ${vcf_file} -dataDir ${projectDir}/../${params.data_directory} -csvStats ${vcf_file.baseName} > ${vcf_file.baseName}.ann.vcf
    """
}

workflow {

    if (params.skip_vcf) {
        def vcf_files = Channel.from(params.vcf_files) 
        snpEff(vcf_files)
    } else {
        bgenToVCF(grouped_files) |
        cleanVCF |
        snpEff
    }
}