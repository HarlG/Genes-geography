## Test initial git commit ##
from pysam import VariantFile
import numpy as np
from sklearn import decomposition
import pandas as pd

vcf_filename = "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
panel_filename = "phase1_integrated_calls.20101123.ALL.panel"


genotypes = []
samples = []
variant_ids = []

#VCF reader
with VariantFile(vcf_filename) as vcf_reader:
    counter = 0
    # iterate through variants in VCF
    for record in vcf_reader:
        counter +=1
        # if statement below only reads every 100 variants to speed up program.
        # This should not disrupt results due to linkage disequilibrium
        # Can be commented out if whole VCF file to be read. Will be option later when args parse functionality added.
        if counter % 100 == 0:
            # Collects allele and sample information
            alleles = [record.samples[x].allele_indices for x in record.samples]
            samples = [sample for sample in record.samples]
            genotypes.append(alleles)
            variant_ids.append(record.id)
        # Print percent completed as program runs
        if counter % 4943 == 0:
            print(f'Completed: {round(100 * counter/494328)}%')
        #if counter >= 10000:
        #    break

# Panel file reader i.e. population codes
with open(panel_filename) as panel_file:
    labels = {} #{sample_id: population code}
    for line in panel_file:
        line = line.strip().split('\t')
        labels[line[0]] = line[1]




genotypes = np.array(genotypes)
print(genotypes.shape)

matrix = np.count_nonzero(genotypes, axis =2)
matrix = matrix.T
print(matrix.shape)

#pca = decomposition.PCA(n_components=2)
#pca.fit(matrix)
#print(pca.singular_values_)
#to_plot = pca.transform(matrix)
#print(to_plot.shape)


# Data put in dataframe. Outputs as CSV
df = pd.DataFrame(matrix, columns = variant_ids, index = samples)
df['Population code'] = df.index.map(labels)
df.to_csv("matrix.csv")