context("single variant tests")
library(SeqVarTools)

test_that("assocTestMM2", {
    svd <- .testData()
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel2(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestMM2(iterator, nullmod, verbose=FALSE)
    restoreFilter(iterator)
    expect_equal(unique(assoc$variant.id), seqGetData(svd, "variant.id"))
    seqClose(svd)
})
