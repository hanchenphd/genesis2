
setGeneric("createDesignMatrix2", function(x, ...) standardGeneric("createDesignMatrix2"))

setMethod("createDesignMatrix2",
          "data.frame",
          function(x, outcome, covars=NULL) {

              if (!is.null(covars)) {
                  model.formula <- as.formula(paste(outcome, "~", paste(covars, collapse="+")))
                  # allow interactions
                  covars <- unique(unlist(strsplit(covars,"[*:]")))
              } else {
                  model.formula <- as.formula(paste(outcome, "~", 1))
              }
              x <- x[, c(outcome, covars), drop=FALSE]
              x <- x[complete.cases(x),,drop=FALSE]
              
              # outcome vector
              y <- x[,outcome]
              # create design matrix    
              X <- model.matrix(model.formula, data=x)
              # check for columns of all the same value (except the intercept)
              dropcol <- append(FALSE, apply(X[,-1,drop=FALSE], 2, var) == 0)
              if (sum(dropcol) > 0) {
                  message("Covariates ",paste(colnames(X)[dropcol], collapse = ", "), " have only 1 value: they have been removed from the model")
                  X <- X[,!dropcol,drop=FALSE]
              }
              list(y=y, X=X)
          })

setMethod("createDesignMatrix2",
          "AnnotatedDataFrame",
          function(x, outcome, covars=NULL, sample.id=NULL) {
              x <- pData(x)
              if (!is.null(sample.id)) {
                  x <- x[x$sample.id %in% sample.id,]
              }
              rownames(x) <- x$sample.id
              createDesignMatrix2(x, outcome, covars)
          })
