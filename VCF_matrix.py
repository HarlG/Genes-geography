## Test initial git commit ##
from pysam import VariantFile


vcf_filename = "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
panel_filename = "phase1_integrated_calls.20101123.ALL.panel"

with VariantFile(vcf_filename) as vcf_reader:
    for record in vcf_reader:
        alleles = [record.samples[x].allele_indices for x in record.samples]
        samples = [sample for sample in record.samples]
        print(alleles)
        print(samples)