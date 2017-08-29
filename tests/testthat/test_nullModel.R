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
