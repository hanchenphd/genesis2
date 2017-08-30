
## account for sex here?
.alleleFreq <- function(geno) {
    0.5*colMeans(geno, na.rm=TRUE)
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

.testData <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    sample.id <- seqGetData(gds, "sample.id")
    df <- data.frame(sample.id=sample.id,
                     sex=sample(c("M","F"), replace=TRUE, length(sample.id)),
                     age=rnorm(length(sample.id), mean=50, sd=10),
                     outcome=rnorm(length(sample.id), mean=10, sd=5),
                     status=rbinom(length(sample.id), size=1, prob=0.4),
                     stringsAsFactors=FALSE)
    sampleData(gds) <- AnnotatedDataFrame(df)
    gds
}
