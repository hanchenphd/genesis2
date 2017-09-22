setGeneric("assocTestMM2", function(gdsobj, ...) standardGeneric("assocTestMM2"))

## arguments to add: test (Wald or Score), ivars
## do we want the ivars.return.betaCov option?
## do we want to make imputing to the mean optional?
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
                  
                  # allele frequency
                  freq <- .alleleFreq(geno)

                  # take note of number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # mean impute missing values
                  if (any(n.obs < nrow(geno))) {
                      geno <- .meanImpute(geno, freq)
                  }

                  # do the test

                  res[[i]] <- cbind(var.info, n.obs, freq)
                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=verbose)
              }

              do.call(rbind, res)
          })
