
# https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/



library(DGEobj.utils)

# Simulate some data
counts <- trunc(matrix(runif(6000, min=0, max=2000), ncol=6))
geneLength <- rowMeans(counts)

# Matrix: N genes x M Samples

# TMM normalized Log2FPKM
Log2FPKM <- convertCounts(counts,
                          unit       = "fpkm",
                          geneLength = geneLength,
                          log        = TRUE,
                          normalize  = "tmm")

# Non-normalized CPM (not logged)
RawCPM <- convertCounts(counts,
                        unit      = "CPM",
                        log       = FALSE,
                        normalize = "none")

# TPM
TPM <- convertCounts(
          counts,
          unit  = "TPM",
          geneLength,
          log = FALSE,
          normalize = "none"
        )


# in case imputation is necessary
library(impute)
df_1_impute <- impute.knn(data.matrix(df_1),k = 10)$data

