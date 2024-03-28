#Last modified on 01/30/2014 

ROSE.eval <- function(response, predictors, data, learner, acc.measure="auc", extr.pred=NULL, method.assess="holdout", K=1, B=100, control.rose=list(), control.learner=list(), control.predict=list(), control.accuracy=list(), trace=FALSE, subset=options("subset")$subset, na_action=options("na.action")$na.action, seed) {
  
  # Check arguments: response, predictors, and learner are mandatory 
  if(missing(response)) 
    stop("Response variable is required.\n")
  if(missing(predictors)) 
    stop("Predictor variables are required.\n")
  if(missing(learner)) 
    stop("Argument 'learner' is missing, with no default. \n")
  
  # Check if provided learner is "standard" in the sense that it has an associated predict method with arguments "object" and "newdata"
  func.name <- as.character(substitute(learner))
  if(any(methods(class=func.name)==paste("predict.",func.name,sep=""))) 
    flg.learner <- 1
  else 
    flg.learner <- 0
  
  mc <- match.call()
  
  # Catch the original data.frame/variables
  data.orig <- sapply(predictors, function(x) get(x, envir=data))
  cn.order.orig <- names(data.orig)
  
  # Keep formula unchanged for the learner
  formula.learn <- predictors
  
  if(missing(data)) {
    lst.model.frame <- list(response=response, predictors=predictors, data=NULL, subset=subset, na_action=na_action)
  } else {
    lst.model.frame <- list(response=response, predictors=predictors, data=data, subset=subset, na_action=na_action)
  }
  
  # Create data set for ROSE and prediction
  mf <- do.call(data.frame, lst.model.frame)
  data.st <- data.frame(mf)
  y <- data.st[, 1]
  
  if(trace) {
    ind <- ifelse(B < 50, 1, ifelse(B < 500, 10, 100))
    cat("Iteration:", "\n")
  }
  
  if(method.assess == "holdout") {
    method.assess <- "BOOT"
    B <- 1
  }
  
  if(method.assess == "BOOT") {
    if(trace) max.ind <- floor(B/ind)*ind
    acc.vec <- numeric(B) 
    
    if(flg.learner) {
      # Functions with "standard" behavior
      for(i in 1:B) {
        data.rose <- do.call(ROSE, c(list(response=response, predictors=predictors, data=data.st), control.rose))$data
        fit <- do.call(learner, c(list(response=response, predictors=predictors, data=data.rose), control.learner))
        pred <- do.call(predict, c(list(object=fit, newdata=data.st), control.predict))
        if(!is.null(extr.pred)) pred <- extr.pred(pred)
        acc.vec[i] <- do.call(fun.accuracy, c(list(response=y, predicted=pred), control.accuracy))[[pos.accuracy]]
        if(trace) if(i %% ind == 0) {if(i != max.ind) cat(i, ", ", sep="") else cat(i, "\n", sep="")}  
      }
    } else {
      # User-defined functions with "non-standard" behavior
      for(i in 1:B) {
        data.rose <- do.call(ROSE, c(list(response=response, predictors=predictors, data=data.st), control.rose))$data
        pred <- do.call(learner, c(list(data=data.rose, newdata=data.st[, -1]), control.learner))
        acc.vec[i] <- do.call(fun.accuracy, c(list(response=y, predicted=pred), control.accuracy))[[pos.accuracy]]
        if(trace) if(i %% ind == 0) {if(i != max.ind) cat(i, ", ", sep="") else cat(i, "\n", sep="")}  
      }
    }
  } else {
    # Leave K out CV
    pred <- y.cp <- numeric(0)
    
    if(trace) max.ind <- floor(B/ind)*ind
    
    if(K %% 1 != 0) stop("Leave K out CV: K must be an integer\n")
    n.g <- K
    
    if(length(data.st[,1]) %% n.g == 0) {
      K <- length(data.st[,1]) / n.g
      ind.g <- sample(rep(1:K, n.g))
    } else {
      K <- floor(length(data.st[,1]) / n.g) + 1
      n.g.remain <- length(data.st[,1]) - floor(length(data.st[,1]) / n.g) * n.g
      message(paste("\nLeave K out CV: the sample size is not a multiple of K. \nThe routine has automatically created", K-1, "subsets of size", n.g, "and one subset of size", n.g.remain,"."))
      ind.g <- sample(c(rep(1:(K-1), n.g), rep(K,n.g.remain) ))
    }
    
    B <- K
    
    if(flg.learner) {
      # Functions with "standard" behavior
      for(i in 1:B) {
        data.rose <- do.call(ROSE, c(list(response=response, predictors=predictors, data=data.st[-which(ind.g==i),]), control.rose))$data
        fit <- do.call(learner, c(list(response=response, predictors=predictors, data=data.rose), control.learner))
        predi <- do.call(predict, c(list(object=fit, newdata=data.st[which(ind.g==i),]), control.predict))
        if(!is.null(extr.pred)) predi <- extr.pred(predi)
        pred <- c(pred, predi)
        if(trace) if(i %% ind == 0) {if(i != max.ind) cat(i, ", ", sep="") else cat(i, "\n", sep="")}  
        y.cp <- c(y.cp, y[which(ind.g==i)])
      }
      acc.vec <- do.call(fun.accuracy, c(list(response=y.cp, predicted=pred), control.accuracy))[[pos.accuracy]]
    } else {
      # User-defined functions with "non-standard" behavior
      for(i in 1:B) {
        data.rose <- do.call(ROSE, c(list(response=response, predictors=predictors, data=data.st[-which(ind.g==i),]), control.rose))$data
        predi <- do.call(learner, c(list(data=data.rose, newdata=data.st[which(ind.g==i), -1]), control.learner))
        pred <- c(pred, predi)
        if(trace) if(i %% ind == 0) {if(i != max.ind) cat(i, ", ", sep="") else cat(i, "\n", sep="")}  
        y.cp <- c(y.cp, y[which(ind.g==i)])
      }
      acc.vec <- do.call(fun.accuracy, c(list(response=y.cp, predicted=pred), control.accuracy))[[pos.accuracy]]
    }
  }
  
  out <- list(Call = mc, method=method.assess, measure = acc.measure, acc = acc.vec)
  class(out) <- "ROSE.eval"
  out
}


##print method for ROSE.eval
print.ROSE.eval <- function(x, ...) 
{
		if (x$method =="BOOT")		method <- "Bootstrap"
		if (x$method =="LKOCV")		method <- "Leave K out cross-validation"
		if (x$method =="holdout")	method <- "Holdout"

	cat("\n")
	cat("Call: \n")
	print(x$Call)
	cat("\n")

		if (method == "Bootstrap")
			cat( paste(method, " estimate of ", x$measure, " on ", length(x$acc), " samples: \n ", sep="") )
		else
			cat( paste(method, " estimate of ", x$measure, ": ", sep="") )

	cat(sprintf("%.3f",x$acc),"\n")
}

###summary method for ROSE.eval
summary.ROSE.eval <- function(object, ...) 
{
	acc<-object$acc
		if (length(acc) > 1) acc <- summary(acc) 
	LST <- list( call=object$Call, method=object$method, measure=object$measure, acc=acc ) 
	class(LST) <- "summary.ROSE.eval"
	LST
}

###print method for summary
print.summary.ROSE.eval <- function(x, ...) 
{
	cat("\n")
	cat("Call: \n")
	print(x$call)
	cat("\n")

		if (x$method =="BOOT") method <- "Bootstrap"
		if (x$method =="LKOCV") method <- "Leave K out cross-validation"
		if (x$method =="holdout") method <- "Holdout"

		if(x$method !="BOOT")
			cat( paste(method, " estimate of ", x$measure, ": ", sprintf("%.3f",x$acc),"\n", sep="") )
		else
		{
			cat( "Summary of bootstrap distribution of auc: \n" )
			print(x$acc)
			cat("\n")
		}
}
