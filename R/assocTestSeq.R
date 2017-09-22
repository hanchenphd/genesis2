setGeneric("assocTestSeq2", function(gdsobj, ...) standardGeneric("assocTestSeq2"))

## should allele frequency filter be an argument, or set ahead of time?
## if an argument, can make sure freq is calcuated on sample set used
## compare seqSetFilterCond to calcluating MAF separately
## previously, returned rows in variantInfo for high-freq variants but set weight to 0
##  - why do this instead of just filtering?

## arguments to add:
## test (Burden or SKAT)
## burden.test (Score or Wald)
## rho
## pval.method (davies, kuonen, liu)

setMethod("assocTestSeq2",
          "SeqVarData",
          function(gdsobj, nullModel,
                   weight.beta=c(1,1), weight.user=NULL,
                   verbose=TRUE) {

              # check argument values

              # results
              res <- list()
              res.var <- list()
              i <- 1
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)
                  
                  geno <- expandedAltDosage(gdsobj)
                  
                  # allele frequency
                  freq <- .alleleFreq(geno)

                  # if we are filtering on freq in function, do it here
                  # then subset var.info, geno, freq

                  # take note of number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)
                  
                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)
                  
                  # mean impute missing values
                  if (any(n.obs < nrow(geno))) {
                      geno <- .meanImpute(geno, freq)
                  }

                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- .weightFromFreq(freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- variantData(gdsobj)[[weight.user]][expandedVariantIndex(gdsobj)]
                      # if we subset based on frequency, will need to do it again here
                  }

                  # do the test

                  res[[i]] <- cbind(n.site, n.alt, n.sample.alt)
                  res.var[[i]] <- cbind(var.info, n.obs, freq, weight)
                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=verbose)
              }

              list(results=do.call(rbind, res), variantInfo=res.var)
          })
