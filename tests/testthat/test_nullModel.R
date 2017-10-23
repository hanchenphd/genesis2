context("null model")
library(Biobase)

test_that("design matrix from data.frame", {
    dat <- data.frame(a=rnorm(10),
                      b=c(rep("a",5), rep("b", 5)),
                      c=sample(1:10, 10, replace=TRUE),
                      d=rep(1, 10))
    dm <- createDesignMatrix2(dat, outcome="a")
    expect_equal(dm$y, dat$a)
    expect_equal(ncol(dm$X), 1)
    expect_true(all(dm$X[,1] == 1))
    dm <- createDesignMatrix2(dat, outcome="a", covars="b")
    expect_equal(colnames(dm$X)[-1], "bb")
    dm <- createDesignMatrix2(dat, outcome="a", covars=c("b", "c", "b:c"))
    expect_equal(colnames(dm$X)[-1], c("bb", "c", "bb:c"))
    expect_message(createDesignMatrix2(dat, outcome="a", covars="d"), "removed from the model")
})

test_that("design matrix from AnnotatedDataFrame", {
    dat <- data.frame(sample.id=sample(letters, 10),
                      a=rnorm(10),
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    dm <- createDesignMatrix2(dat, outcome="a", covars="b")
    expect_equal(dm$y, dat$a)
    expect_equal(rownames(dm$X), dat$sample.id)
    keep <- dat$sample.id[c(TRUE,FALSE)]
    dm <- createDesignMatrix2(dat, outcome="a", covars="b", sample.id=keep)
    expect_equal(dm$y, dat$a[c(TRUE,FALSE)])
    expect_equal(rownames(dm$X), keep)
})

test_that("null model", {
    dat <- data.frame(sample.id=sample(letters, 10),
                      a=rnorm(10),
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    keep <- dat$sample.id[c(TRUE,FALSE)]
    nm <- fitNullModel2(dat, outcome="a", covars="b", sample.id=keep, verbose=FALSE)
    expect_equal(rownames(nm$model.matrix), keep)
    expect_equal(nm$workingY, dat$a[c(TRUE,FALSE)])
})

test_that("index list", {
    x <- rep(letters[1:3], each=3)
    expect_equal(list(a=1:3, b=4:6, c=7:9), .indexList(x))
    expect_equal(list(a=1:3), .indexList(rep("a", 3)))
})



test_that("fitNullMM2 matches fitNulMM", {
    library(SeqVarTools)
    library(Biobase)
    svd <- .testData()
    grm <- .testGRM(svd)
    lmm.genesis <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    nullmod <- fitNullModel2(svd, outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)

    expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-9))
    expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-9))
    expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-9))
    expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-9))
    expect_true(all(abs(nullmod$logLikR - lmm.genesis$logLikR) < 1e-9))
    expect_true(all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-9))
    expect_true(all(nullmod$workingY == lmm.genesis$workingY))
    expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
    expect_true(all(abs(nullmod$varComp - lmm.genesis$varComp) < 1e-9))
    expect_true(all(abs(nullmod$varCompCov - lmm.genesis$varCompCov) < 1e-9)) 
    expect_equal(nullmod$family$family, lmm.genesis$family$family)
    expect_true(all(nullmod$zeroFLAG == lmm.genesis$zeroFLAG))
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-9))
    expect_true(all(abs(nullmod$RSS - lmm.genesis$RSS) < 1e-9))

    seqClose(svd)
})

test_that("fitNullMM2 matches fitNulMM - group", {
    library(SeqVarTools)
    library(Biobase)
    svd <- .testData()
    grm <- .testGRM(svd)
    lmm.genesis <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, group.var="status", verbose=FALSE)
    nullmod <- fitNullModel2(svd, outcome="outcome", covars=c("sex", "age"), covMatList=grm, group.var="status", verbose=FALSE)

    expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-9))
    expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-9))
    expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-9))
    expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-9))
    expect_true(all(abs(nullmod$logLikR - lmm.genesis$logLikR) < 1e-9))
    expect_true(all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-9))
    expect_true(all(nullmod$workingY == lmm.genesis$workingY))
    expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
    expect_true(all(abs(nullmod$varComp - lmm.genesis$varComp) < 1e-9))
    expect_equal(nullmod$varCompCov, lmm.genesis$varCompCov) 
    expect_equal(nullmod$family$family, lmm.genesis$family$family)
    expect_true(all(nullmod$zeroFLAG == lmm.genesis$zeroFLAG))
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-9))
    expect_true(all(abs(nullmod$RSS - lmm.genesis$RSS) < 1e-9))

    seqClose(svd)
})

test_that("fitNullMM2 matches fitNulMM - binary", {
    library(SeqVarTools)
    library(Biobase)
    svd <- .testData()
    grm <- .testGRM(svd)
    lmm.genesis <- GENESIS::fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
    nullmod <- fitNullModel2(svd, outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)

    expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-6))
    expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-6))
    expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-6))
    expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-6))
    expect_true(all(abs(nullmod$logLikR - lmm.genesis$logLikR) < 1e-6))
    expect_true(all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-6))
    expect_true(all(abs(nullmod$workingY - lmm.genesis$workingY) < 1e-6))
    expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
    expect_true(all(abs(nullmod$varComp - lmm.genesis$varComp) < 1e-6))
    expect_true(all(abs(nullmod$varCompCov - lmm.genesis$varCompCov) < 1e-6)) 
    expect_equal(nullmod$family$family, lmm.genesis$family$family)
    expect_true(all(nullmod$zeroFLAG == lmm.genesis$zeroFLAG))
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-6))
    #expect_true(all(abs(nullmod$RSS - lmm.genesis$RSS) < 1e-6))

    seqClose(svd)
})

test_that("fitNullMM2 matches fitNullReg", {
    library(SeqVarTools)
    library(Biobase)
    svd <- .testData()
    lmm.genesis <- GENESIS::fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    nullmod <- fitNullModel2(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)

    expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-9))
    expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-9))
    expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-9))
    expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-9))
    expect_true(all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-9))
    expect_true(all(nullmod$workingY == lmm.genesis$workingY))
    expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
    expect_equal(nullmod$family$family, lmm.genesis$family$family)

    seqClose(svd)
})

test_that("fitNullMM2 matches fitNullReg - binary", {
    library(SeqVarTools)
    library(Biobase)
    svd <- .testData()
    lmm.genesis <- GENESIS::fitNullReg(sampleData(svd), outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
    nullmod <- fitNullModel2(svd, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)

    expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-9))
    expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-9))
    expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-9))
    expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-9))
    expect_true(all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-9))
    expect_true(all(nullmod$workingY == lmm.genesis$workingY))
    expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
    expect_equal(nullmod$family$family, lmm.genesis$family$family)

    seqClose(svd)
})
