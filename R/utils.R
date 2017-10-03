
.alleleFreq <- function(gdsobj, geno) {

    # check sex
    sex <- validateSex(gdsobj)
    if (is.null(sex)) {
        return(0.5*colMeans(geno, na.rm=TRUE))
    }

    # check chromosome
    chr <- chromWithPAR(gdsobj)
    X <- chr %in% "X"
    Y <- chr %in% "Y"
    auto <- !X & !Y

    # allele frequency vector
    freq <- rep(NA, ncol(geno))

    # autosomes
    if (any(auto)) {
        freq[auto] <- 0.5*colMeans(geno[, auto, drop=FALSE], na.rm=TRUE)
    }

    # X chrom
    if (any(X)) {
        female <- sex %in% "F"
        male <- sex %in% "M"
        F.count <- colSums(geno[female, X, drop=FALSE], na.rm=TRUE)
        F.nsamp <- colSums(!is.na(geno[female, X, drop=FALSE]))
        M.count <- 0.5*colSums(geno[male, X, drop=FALSE], na.rm=TRUE)
        M.nsamp <- colSums(!is.na(geno[male, X, drop=FALSE]))
        freq[X] <- (F.count + M.count)/(2*F.nsamp + M.nsamp)
    }

    # Y chrom
    if (any(Y)) {
        male <- sex %in% "M"
        freq[Y] <- 0.5*colMeans(geno[male, Y, drop=FALSE], na.rm=TRUE)
    }

    freq
}

.meanImpute <- function(geno, freq) {
    miss.idx <- which(is.na(geno))
    miss.var.idx <- ceiling(miss.idx/nrow(geno))
    geno[miss.idx] <- 2*freq[miss.var.idx]
    geno
}

.weightFromFreq <- function(freq, weight.beta) {
    freq <- pmin(freq, 1-freq)
    dbeta(freq, weight.beta[1], weight.beta[2])
}
