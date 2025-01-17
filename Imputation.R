library(Rmagic)
library(trqwe)

impute <- function(Sample) {
    count_matrix = mcreadRDS(paste0("../Data_generated/",Sample,"_TPM.rds"), mc.cores=4)
    count_matrix = t(count_matrix)
    gc()
    count_matrix = magic(count_matrix, n.jobs = -2, solver='exact')
    gc()
    mcsaveRDS(t(count_matrix$result), paste0("../Data_generated/",Sample,"_imputed_exact.rds"), mc.cores = 4)
}