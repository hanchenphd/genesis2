setGeneric("assocTestMM2", function(gdsobj, ...) standardGeneric("assocTestMM2"))

setMethod("assocTestMM2",
          "SeqVarData",
          function(gdsobj, nullModel, verbose=TRUE) {
              # results
              res <- list()
              i <- 1
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)
                  
                  geno <- expandedAltDosage(gdsobj)
                  
                  # take note of number of non-missing samples
                  n.miss <- colSums(is.na(geno))
                  n.obs <- nrow(geno) - n.miss
                  
                  # allele frequency (account for sex?)
                  freq <- 0.5*colMeans(geno, na.rm=TRUE)

                  # mean impute missing values
                  if (sum(n.miss) > 0) {
                      geno <- .meanImpute(geno, freq)
                  }

                  # do the test

                  res[[i]] <- cbind(var.info, n.obs, freq)
                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=verbose)
              }

              do.call(rbind, res)
          })

.meanImpute <- function(geno, freq) {
    miss.idx <- which(is.na(geno))
    miss.var.idx <- ceiling(miss.idx/nrow(geno))
    geno[miss.idx] <- 2*freq[miss.var.idx]
    geno
}
