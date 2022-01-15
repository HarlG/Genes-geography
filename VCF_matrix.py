## Test initial git commit ##
from pysam import VariantFile
import numpy as np


vcf_filename = "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
panel_filename = "phase1_integrated_calls.20101123.ALL.panel"


genotypes = []
samples = []

with VariantFile(vcf_filename) as vcf_reader:
    counter = 0
    for record in vcf_reader:
        alleles = [record.samples[x].allele_indices for x in record.samples]
        samples = [sample for sample in record.samples]
        genotypes.append(alleles)
        counter += 1
        if counter >= 100:
            break
        
genotypes = np.array(genotypes)
print(genotypes.shape)
