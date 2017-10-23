
setGeneric("fitNullModel2", function(x, ...) standardGeneric("fitNullModel2"))

setMethod("fitNullModel2",
          "AnnotatedDataFrame",
          function(x, outcome,
                   covars = NULL,
                   covMatList = NULL,
                   group.var = NULL,
                   sample.id = NULL,
                   family = "gaussian",
                   start = NULL,
                   AIREML.tol = 1e-6,
                   maxIter = 100,
                   dropZeros = TRUE,
                   verbose = TRUE) {
              desmat <- createDesignMatrix2(x, outcome, covars, sample.id)
              group.idx <- .indexList(x[[group.var]])
              fitNullModel(y=desmat$y, X=desmat$X, covMatList=covMatList, group.idx=group.idx,
                           family=family, start=start, AIREML.tol=AIREML.tol, maxIter=maxIter,
                           dropZeros=dropZeros, verbose=verbose)
          })

setMethod("fitNullModel2",
          "SeqVarData",
          function(x, ...) {
              fitNullModel2(sampleData(x), ...)
          })

.indexList <- function(x) {
    groups <- unique(x)
    lapply(setNames(groups, groups), function(g) which(x == g))
}
