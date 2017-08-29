context("single variant tests")
library(SeqVarTools)
library(Biobase)

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

test_that("assocTestMM2", {
    svd <- .testData()
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel2(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestMM2(iterator, nullmod, verbose=FALSE)
    seqClose(svd)
})
