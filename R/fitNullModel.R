
setGeneric("fitNullModel2", function(x, ...) standardGeneric("fitNullModel2"))

setMethod("fitNullModel2",
          "AnnotatedDataFrame",
          function(x, outcome, covars=NULL, sample.id=NULL, verbose=TRUE, ...) {
              desmat <- createDesignMatrix2(x, outcome, covars, sample.id)
              fitNullModel(y=desmat$y, X=desmat$X, verbose=verbose, ...)
          })

setMethod("fitNullModel2",
          "SeqVarData",
          function(x, ...) {
              fitNullModel2(sampleData(x), ...)
          })
