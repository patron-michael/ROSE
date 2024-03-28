#Last modified on 01/30/2014 
######################################################################
.onAttach <- function(libname,pkgname){
   packageStartupMessage("Loaded ROSE ", as.character(packageDescription("ROSE")[["Version"]]),"\n")
}
######################################################################

######################################################################
#new.ovun.sample main function
######################################################################
new.ovun.sample <- function(response_var, predictor_vars, data, method="both", N, p=0.5, subset=options("subset")$subset, na.action=options("na.action")$na.action, seed)
{

	###checks
		if(is.null(predictor_vars)) 
			stop("predictors are reaquired.\n")
  
    if(is.null(response_var)) 
      stop("variables are reaquired.\n")

	method <- match.arg(method, choices=c("both", "under", "over"))
		if( !method%in%c("both", "over", "under") ) 
			stop("Method must be 'both', 'over', or 'under'.\n")
	###
	Call <- match.call()
	m <- match(c("response_var", "predictor_vars", "data","method","N", "p", "seed", "subset", "na.action"), names(Call), 0L)
	Call1 <- Call[c(1L, m)]
	Call1[[1L]] <- new.omnibus.balancing
	res <- eval(Call1)
	out <- list(Call=match.call(), method=method, data=res$data)
	class(out) <- "new.ovun.sample"
	out
}

##print method for new.ovun.sample
print.new.ovun.sample <- function(x, ...) 
{
	cat("\n")
	cat("Call: \n")
	print(x$Call)
	Method <- switch(match.arg(x$method, choices=c("both", "under", "over")),
							both="combination of over- and under-sampling",
							under="undersampling",
							over="oversampling"
						 )
	cat("\n")
	cat("Data balanced by", Method,"\n")
	cat("\n")
	print(x$data)
}

###summary method for new.ovun.sample
summary.new.ovun.sample <- function(object, ...) 
{
	out <- list( Call=object$Call, Summary=summary(object$data), method=object$method )
	class(out) <- "summary.new.ovun.sample"
	out
}

###print method for summary new.ovun.sample
print.summary.new.ovun.sample <- function(x, ...) 
{
	cat("\n")
	cat("Call: \n")
	print(x$Call)
	cat("\n")

	Method <- switch(match.arg(x$method, choices=c("both", "under", "over")),
							both="combination of over- and under-sampling",
							under="undersampling",
							over="oversampling"
						 )

	cat("Summary of data balanced by", Method ,"\n")
	cat("\n")
	print(x$Summary)
}

######################################################################
##function that provides a formula with non tranformed variables only
######################################################################
##this function is NOT exported
adj.formula <- function(formula, data)
{
	if( missing(data) )
		frml.env <- environment(formula)
	else
		frml.env <- data

	formula <- terms(formula, data = frml.env)
	vars <- attr(formula, "variables")
	vars <- sapply(vars, function(x) paste(deparse(x,width.cutoff=500), collapse=' '))[-1L]
	#remove all characters before either ( or /
	vars <- sub("*.*[(/]","", vars)
	#remove all characters after either ^ or )
	vars <- sub("['^')].*","", vars)
	vars <- unique(vars)
	formula <- as.formula(paste(vars[1], "~", paste(vars[-1], collapse= "+")))
	attr(formula, "variables") <- vars
	formula
}


######################################################################
#This function is the wrapper for all the implemented data balancing remedies 
######################################################################
##this function is NOT exported
new.omnibus.balancing <- function(response_var, predictor_vars, data, method, subset, na.action, N, p=0.5, seed, hmult.majo=1, hmult.mino=1)
{

  ### Checks and argument parsing
  if(is.null(predictor_vars)) 
    stop("Predictor variables are required.\n")
  if(is.null(response_var)) 
    stop("Response variable is required.\n")
  if(missing(method))
    method <- "both"
  if((method == "under" || method == "over") && !missing(N) && !missing(p))
    stop("Too many arguments. Need to specify either N or p.\n")
  
  # response_var <- as.character(substitute(response_var))
  # predictor_vars <- as.character(substitute(predictor_vars))
  
  response <- data[[response_var]]
  predictors <- data[, predictor_vars, drop = FALSE]
  
  # Handle missing arguments
  if(missing(subset))
    subset <- options("subset")$subset
  if(missing(na.action))
    na_action <- options("na.action")$na.action
  
  # Extracting necessary variables and performing necessary checks
  n <- length(response)
  d <- NCOL(predictors)
  
  classy <- class(response)
  response <- factor(response)
  T <- table(response)
  classx <- sapply(as.data.frame(predictors), class)
  
  if(n < 2) 
    stop("Too few observations.\n")  
  
  if(length(T) > 2)
    stop("The response variable must have 2 levels.\n")
  else if(length(T) == 1)
    stop("The response variable has only one class.\n")
  
  if(p < 0 || p > 1) 
    stop("p must be in the interval 0-1.\n")
  
  majoY <- levels(response)[which.max(T)]
  minoY <- levels(response)[which.min(T)]
  
  ind.mino <- which(response == minoY)
  ind.majo <- which(response == majoY)
  
  if(!missing(seed)) 
    set.seed(seed)
  
  # Handling the selected method
  data.obj <- switch(method,
                     both = ou.sampl(n, N, p, ind.majo, majoY, ind.mino, minoY, classy, predictors),
                     over = over.sampl(n, N, p, ind.majo, ind.mino, majoY, minoY, response, classy, predictors),
                     under = under.sampl(n, N, p, ind.majo, majoY, ind.mino, minoY, response, classy, predictors),
                     rose = rose.sampl(n, N, p, ind.majo, majoY, ind.mino, minoY, response, classy, predictors, classx, d, T, hmult.majo, hmult.mino)
  )
  
  data.out <- data.obj$data.out
  ynew <- data.obj$ynew
  Xnew <- data.obj$Xnew
  
  # Re-positioning columns if necessary
  if(!missing(data)) {
    colnames(data.out) <- colnames(data)[colnames(data) %in% colnames(data.out)]
    data.out <- cbind(ynew, Xnew)
  }
  
  # TODO: adapt for have the truth names
  names(data.out) <- names(data)
  list(data = data.out, call = match.call())
}

######################################################################
#Combination of over and under sampling
######################################################################
##this function is NOT exported
ou.sampl <- function(n, N, p, ind.majo, majoY, ind.mino, minoY, classy, X)
{

		if( missing(N) )
			N <- n
	#number of new minority class examples
	n.mino.new <- sum(rbinom(N, 1, p))
	#number of new majority class examples
	n.majo.new <- N-n.mino.new

	id.majo.new <- sample(ind.majo, n.majo.new, replace=TRUE)
	id.mino.new <- sample(ind.mino, n.mino.new, replace=TRUE)

	#create X
	Xnew <- data.frame(X[c(id.majo.new, id.mino.new),])
	#create  y
		if( classy%in%c("character", "integer", "numeric") )
			ynew <- as.vector( c(rep(majoY, n.majo.new), rep(minoY, n.mino.new)), mode=classy )
		if( classy=="factor" )  
			ynew <- factor( c(rep(majoY, n.majo.new), rep(minoY, n.mino.new)), levels=c(majoY, minoY) )

	data.out <- data.frame(ynew, Xnew)
	rownames(data.out) <- 1:N

	list(data.out=data.out, ynew=ynew, Xnew=Xnew)
}

######################################################################
#Under sampling
######################################################################
##this function is NOT exported
under.sampl <- function(n, N, p, ind.majo, majoY, ind.mino, minoY, y, classy, X)
{

	n.mino.new <- sum(y == minoY)

		if( missing(N) )
		{
				# Determination of N and n.majo in version 0.0.2
				if( p<n.mino.new/n ) 
					warning("non-sensible to specify p smaller than the actual proportion of minority class examples in the original sample.\n")
			#theoretical n.majo
			n.majo <- round( (1-p)*n.mino.new/p )
			#estimated n.majo
			n.majo.new <- sum( rbinom(n.mino.new+n.majo, 1, 1-p) )
			#final sample size
			N <- n.majo.new + n.mino.new
		}
		else
		{
			if(N<n.mino.new)
				stop("N must be greater or equal than the number of minority class examples.\n")
			else
				n.majo.new <- N-n.mino.new
		}

	id.mino.new <- ind.mino
	id.majo.new <- sample(ind.majo, n.majo.new, replace=FALSE)

	#create X
	Xnew <- data.frame(X[c(id.majo.new, id.mino.new),])
	#create  y
		if( classy%in%c("character", "integer", "numeric") )
			ynew <- as.vector( c(rep(majoY, n.majo.new), rep(minoY, n.mino.new)), mode=classy )
		if( classy=="factor" )  
			ynew <- factor( c(rep(majoY, n.majo.new), rep(minoY, n.mino.new)), levels=c(majoY, minoY) )

	data.out <- data.frame(ynew, Xnew)
	rownames(data.out) <- 1:N

	list(data.out=data.out, ynew=ynew, Xnew=Xnew)

}

######################################################################
#Over sampling
######################################################################
over.sampl <- function(n, N, p, ind.majo, ind.mino, majoY, minoY, y, classy, X)
{

	n.majo <- n.majo.new <- sum(y == majoY)
	n.mino <- n-n.majo

		if( missing(N) )
		{
				if( p<n.mino/n ) 
					warning("non-sensible to specify p smaller than the actual proportion of minority class examples in the original sample.\n")
				#theoretical n.mino
				n.mino <- round( p*n.majo/(1-p) )
				#estimated n.mino
				n.mino.new <- sum( rbinom(n.mino+n.majo, 1, p) )
				#final sample size
				N <- n.majo + n.mino.new
		}
		else
		{
				if(N<n)
					stop("N must be greater or equal than the actual sample size.\n")
				else
					n.mino.new <- N-n.majo
		}

	id.majo.new <- ind.majo
	id.mino.new <- sample(ind.mino, n.mino.new, replace=TRUE)

	#create X
	Xnew <- data.frame(X[c(id.majo.new, id.mino.new),])
	#create  y
		if( classy%in%c("character", "integer", "numeric") )
			ynew <- as.vector( c(rep(majoY, n.majo.new), rep(minoY, n.mino.new)), mode=classy )
		if( classy=="factor" )  
			ynew <- factor( c(rep(majoY, n.majo.new), rep(minoY, n.mino.new)), levels=c(majoY, minoY) )

	data.out <- data.frame(ynew, Xnew)
	rownames(data.out) <- 1:N

	list(data.out=data.out, ynew=ynew, Xnew=Xnew)

}

######################################################################
#Rose generation
######################################################################
rose.sampl <- function(n, N, p, ind.majo, majoY, ind.mino, minoY, y, classy, X, classx, d, T, hmult.majo, hmult.mino)
{
		# variables: must be numeric, integer or factor
		if( any( is.na( pmatch(classx, c( "numeric","integer","factor"), duplicates.ok = TRUE ) ) ) ) 
			stop("The current implementation of ROSE handles only continuous and categorical variables.\n")

		if( any(T < 2) ) 
			stop("ROSE needs at least two majority and two minority class examples.\n")

		if( missing(N) )
			N <- n
	#number of new minority class examples
	n.mino.new <- sum(rbinom(N, 1, p))
	#number of new majority class examples
	n.majo.new <- N-n.mino.new

	id.majo.new <- sample(ind.majo, n.majo.new, replace=TRUE)
	id.mino.new <- sample(ind.mino, n.mino.new, replace=TRUE)


	id.num  <- which(classx=="numeric" | classx=="integer")
	d.num   <- d-length( which(classx=="factor") )

	#create  X
	Xnew <- data.frame(X[c(id.majo.new, id.mino.new),])
		if(d.num > 0)  
		{
			 Xnew[1:n.majo.new, id.num] <- rose.real(X[,id.num], hmult=hmult.majo, n=length(ind.majo), q=d.num, ids.class=ind.majo, ids.generation=id.majo.new)
			 Xnew[(n.majo.new+1):N, id.num] <- rose.real(X[,id.num], hmult=hmult.mino, n=length(ind.mino), q=d.num, ids.class=ind.mino, ids.generation=id.mino.new)
		}

	#create  y
		if( classy%in%c("character", "integer", "numeric") )
			ynew <- as.vector( c(rep(majoY, n.majo.new), rep(minoY, n.mino.new)), mode=classy )
		if( classy=="factor" )  
			ynew <- factor( c(rep(majoY, n.majo.new), rep(minoY, n.mino.new)), levels=c(majoY, minoY) )

	data.out <- data.frame(ynew, Xnew)
	rownames(data.out) <- 1:N

	list(data.out=data.out, ynew=ynew, Xnew=Xnew)
}


######################################################################
#function to generate synthetic real data
######################################################################
##This function is NOT exported
rose.real <- function(X, hmult=1, n, q = NCOL(X), ids.class, ids.generation)
{
	X <- data.matrix(X)
	n.new <- length(ids.generation)
	cons.kernel <- (4/((q+2)*n))^(1/(q+4))

		if(q!=1)
			H <- hmult*cons.kernel*diag(apply(X[ids.class,], 2, sd), q)
		else
			H <- hmult*cons.kernel*sd(X[ids.class,])

	Xnew.num <- matrix(rnorm(n.new*q), n.new, q)%*%H
	Xnew.num <- data.matrix(Xnew.num + X[ids.generation,])
	Xnew.num
}

######################################################################
#Wrapper for ROSE
######################################################################
ROSE <- function(response_var, predictor_vars, data, N, p=0.5, hmult.majo=1, hmult.mino=1, subset=options("subset")$subset, na.action=options("na.action")$na.action, seed)
{
  mc <- match.call()
  obj <- new.omnibus.balancing(response_var, predictor_vars, data, subset, na.action, N, p, method="rose", seed, hmult.majo, hmult.mino)
  # out <- list(Call=mc, method="ROSE", data=obj$data)
  # class(out) <- "ROSE"
  # out
}


##print method for ROSE
print.ROSE <- function(x, ...) 
{
	cat("\n")
	cat("Call: \n")
	print(x$Call)
	cat("\n")
	cat("Data balanced by", x$method,"\n")
	cat("\n")
	print(x$data)
}

###summary method for ROSE
summary.ROSE <- function(object, ...) 
{
	out <- list( Call=object$Call, Summary=summary(object$data) )
	class(out) <- "summary.ROSE"
	out
}

###print method for summary ROSE
print.summary.ROSE <- function(x, ...) 
{
	cat("\n")
	cat("Call: \n")
	print(x$Call)
	cat("\n")

	cat("Summary of data balanced by ROSE","\n")
	cat("\n")
	print(x$Summary)
}

