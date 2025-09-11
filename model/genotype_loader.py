import glow
from pyspark.sql import SparkSession


spark = (
    SparkSession.builder
    .appName("VCF Loader")
    # Tell Spark to pull in the Glow JAR
    .config("spark.jars.packages", "io.projectglow:glow-spark3_2.12:2.0.0")
    # Required so Glow can read/write BGZF-compressed files (VCF bgz/gz)
    .config("spark.hadoop.io.compression.codecs", "io.projectglow.sql.util.BGZFCodec")
    .getOrCreate()
)

# Register Glow functions with Spark
glow.register(spark)

vcf_file = "PGS001990/ukb_imp_bct_vars_merged_effects.vcf"

# Example: load VCF
df = (
    spark.read.format("vcf")
    .option("includeSampleIds", True)
    .load(vcf_file)
)


#vcf_reader = vcfpy.Reader.from_path(vcf)

#samples = vcf_reader.header.samples.names
## for testing
#samples = samples[0:5]
#print(samples)

#sample_to_genotypes = { samp: {} for samp in samples }

#print(sample_to_genotypes)

#i = 0
#for record in vcf_reader:
        ## Determine an identifier for the variant
        ## If record.ID is set and not “.”, use that; else make one
    #varid = record.ID
    #print(varid)
    #for call in record.calls[0:5]:
        #samp = call.sample
        ## call.data is an OrderedDict containing FORMAT fields
        ## The genotype is under key “GT”
        #gt = call.data.get('GT')
        #print(samp, gt)

        #if gt == '0/0':
            #gt_score = 0
        #elif gt == '0/1':
            #gt_score = 1
        #elif gt == '1/1':
            #gt_score = 2
        #else:
            #print(f'{gt} ??')

        #sample_to_genotypes[samp][varid[0]] = gt_score

    #if i == 4:
        #print(sample_to_genotypes)
        #exit()
    #i += 1