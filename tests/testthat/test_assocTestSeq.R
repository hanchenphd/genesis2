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


test_that("assocTestSeq2 matches assocTestSeqWindow - Burden, Wald", {
    svd <- .testData()
    nullmod <- GENESIS::fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Wald", verbose=FALSE)
    
    nullmod <- fitNullModel2(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestSeq2(iterator, nullmod, test="Burden", burden.test="Wald", verbose=FALSE)

    var1 <- assoc1$variantInfo[!(assoc1$variantInfo$freq %in% c(0,1)),]
    var2 <- do.call(rbind, assoc2$variantInfo)
    var2 <- var2[!duplicated(paste(var2$variant.id, var2$allele.index)),]
    expect_equal(nrow(var1), nrow(var2))  
    expect_equal(var1$variantID, var2$variant.id)
    expect_equal(var1$allele, var2$allele.index)
    expect_equal(var1$chr, var2$chromosome)
    expect_equal(var1$pos, var2$position)
    expect_equal(var1$n.obs, var2$n.obs)
    expect_equal(var1$freq, var2$freq)
    expect_equal(var1$weight, var2$weight)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    chk <- function(x) {
        paste(x$n.site, format(x$Est, digits=4, scientific=FALSE, trim=TRUE), sep="_")
    }
    expect_true(setequal(chk(res1), chk(res2)))
    
    seqClose(svd)
})

test_that("assocTestSeq2 matches assocTestSeqWindow - Burden, Score", {
    svd <- .testData()
    nullmod <- GENESIS::fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE)
    
    nullmod <- fitNullModel2(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestSeq2(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    chk <- function(x) {
        paste(x$n.site, format(x$Score, digits=4, scientific=FALSE, trim=TRUE), sep="_")
    }
    expect_true(setequal(chk(res1), chk(res2)))
    
    seqClose(svd)
})

test_that("assocTestSeq2 matches assocTestSeqWindow - SKAT", {
    svd <- .testData()
    nullmod <- GENESIS::fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE)
    
    nullmod <- fitNullModel2(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestSeq2(iterator, nullmod, test="SKAT", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    chk <- function(x) {
        fmt <- function(y) format(y, digits=4, scientific=FALSE, trim=TRUE)
        paste(x$n.site, fmt(x$Q_0), fmt(x$pval_0), fmt(x$err_0), sep="_")
    }
    expect_true(setequal(chk(res1), chk(res2)))
    
    seqClose(svd)
})
