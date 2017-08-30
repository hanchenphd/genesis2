context("rare variant tests")
library(SeqVarTools)
library(GenomicRanges)
library(Biobase)

test_that("window", {
    svd <- .testData()
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
    nullmod <- fitNullModel2(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestSeq2(iterator, nullmod, verbose=FALSE)
    nwin <- length(variantRanges(iterator))
    expect_equal(nrow(assoc$results), nwin)
    expect_equal(length(assoc$variantInfo), nwin)
    seqClose(svd)
})

test_that("ranges", {
    svd <- .testData()
    gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
    iterator <- SeqVarRangeIterator(svd, variantRanges=gr, verbose=FALSE)
    nullmod <- fitNullModel2(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestSeq2(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(assoc$results), length(gr))
    expect_equal(length(assoc$variantInfo), length(gr))
    seqClose(svd)
})

test_that("list", {
    svd <- .testData()
    gr <- GRangesList(
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(1e6, 3e6), width=1e6)),
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(3e6, 34e6), width=1e6)))
    iterator <- SeqVarListIterator(svd, variantRanges=gr, verbose=FALSE)
    assoc <- assocTestSeq2(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(assoc$results), length(gr))
    expect_equal(length(assoc$variantInfo), length(gr))
    seqClose(svd)
})

test_that("user weights", {
    svd <- .testData()
    variant.id <- seqGetData(svd, "variant.id")
    weights <- sample(1:10, length(variant.id), replace=TRUE)
    variantData(svd) <- AnnotatedDataFrame(data.frame(variant.id, weights))
    seqSetFilterChrom(svd, include=21, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
    assoc <- assocTestSeq2(iterator, nullmod, weight.user="weights", verbose=FALSE)
    tmp <- do.call(rbind, assoc$variantInfo)[,c("variant.id", "weight")]
    tmp <- tmp[!duplicated(tmp$variant.id),]
    restoreFilter(svd)
    expect_equal(tmp$weight, variantData(svd)$weight)
    seqClose(svd)
})
