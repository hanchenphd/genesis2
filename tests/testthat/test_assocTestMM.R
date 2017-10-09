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


test_that("assocTestMM2 matches assocTestMM - Wald", {
    svd <- .testData()
    grm <- .testGRM(svd)
    seqResetFilter(svd, verbose=FALSE)
    nullmod <- fitNullModel2(svd, outcome="outcome", covMatList=grm, verbose=FALSE)
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    assoc1 <- GENESIS::assocTestMM(svd, nullmod, snp.include=which(snv), verbose=FALSE)
    
    seqSetFilter(svd, variant.sel=snv, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    assoc2 <- assocTestMM2(iterator, nullmod, verbose=FALSE)
    not1 <- setdiff(assoc2$variant.id, assoc1$snpID)
    expect_true(all(assoc2$freq[not1] %in% c(0,1)))
    expect_true(all(is.na(assoc2$Est[not1])))
    assoc2 <- assoc2[!(assoc2$variant.id %in% not1),]
    expect_equal(nrow(assoc1), nrow(assoc2))
    expect_equal(assoc1$snpID, assoc2$variant.id)
    expect_equal(assoc1$chr, assoc2$chromosome)
    expect_equal(assoc1$MAF, pmin(assoc2$freq, 1-assoc2$freq))
    expect_equal(assoc1$Est, assoc2$Est)
    expect_equal(assoc1$SE, assoc2$Est.SE)
    expect_equal(assoc1$Wald.Stat, (assoc2$Wald.Stat)^2)
    expect_equal(assoc1$Wald.pval, assoc2$Wald.pval)
    
    seqClose(svd)
})


test_that("assocTestMM2 matches assocTestMM - Score", {
    svd <- .testData()
    grm <- .testGRM(svd)
    seqResetFilter(svd, verbose=FALSE)
    nullmod <- fitNullModel2(svd, outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    assoc1 <- GENESIS::assocTestMM(svd, nullmod, test="Score", snp.include=which(snv), verbose=FALSE)
    
    seqSetFilter(svd, variant.sel=snv, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    assoc2 <- assocTestMM2(iterator, nullmod, test="Score", verbose=FALSE)
    not1 <- setdiff(assoc2$variant.id, assoc1$snpID)
    assoc2 <- assoc2[!(assoc2$variant.id %in% not1),]
    expect_equal(nrow(assoc1), nrow(assoc2))
    expect_equal(assoc1$Score, assoc2$Score)
    expect_equal(assoc1$Var, (assoc2$Score.SE)^2)
    expect_equal(assoc1$Score.Stat, (assoc2$Score.Stat)^2)
    expect_equal(assoc1$Score.pval, assoc2$Score.pval)
    
    seqClose(svd)
})
