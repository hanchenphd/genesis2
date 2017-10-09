setGeneric("assocTestMM2", function(gdsobj, ...) standardGeneric("assocTestMM2"))

## do we want the ivars.return.betaCov option?
## do we want to make imputing to the mean optional?
setMethod("assocTestMM2",
          "SeqVarData",
          function(gdsobj, nullModel, test = c("Wald", "Score"), ivars = NULL, verbose=TRUE) {
              test <- match.arg(test)
              
              # results
              res <- list()
              i <- 1
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)
                  
                  geno <- expandedAltDosage(gdsobj)
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno)

                  # take note of number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # mean impute missing values
                  if (any(n.obs < nrow(geno))) {
                      geno <- .meanImpute(geno, freq)
                  }

                  # do the test
                  nullPrep <- nullModelTestPrep(nullModel)
                  assoc <- testGenoSingleVar(nullPrep, G=geno, E=ivars, test=test)
                  # set monomorphs to NA - do we want to skip testing these to save time?
                  assoc[freq %in% c(0,1),] <- NA

                  res[[i]] <- cbind(var.info, n.obs, freq, assoc)
                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=verbose)
              }

              do.call(rbind, res)
          })
