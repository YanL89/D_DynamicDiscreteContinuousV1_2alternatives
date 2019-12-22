#' create a design matrix for a utility based discrete model
#' 
#' @description \code{create_X} creates a design matrix X such that X beta is a vector of all the utilities 
#' 
#' @param generic a list of vectors that specify variables who share a generic coefficient for all alternatives
#' @param specific a list of vectors of all parameters with a specific (unique) coefficient for each alternative.
#' @param D the data set
#' @param ASC ASC specification, TRUE, FALSE or vector of id that need ASC. NULL is cast to FALSE
#' @return the design matrix X
#' @export 
create_X = function(generic, specific, D, ASC = NULL){
	# check for ASC
	if(is.null(ASC))
		ASC = FALSE
	specific = augmentSpecific(specific, ASC)
	
  nalt = length(specific)
  n = nrow(D)
  nrow = nalt * n
  npar = length(generic)
  for(i in specific)
    npar = npar + length(i)
  
  X = matrix(0,nrow = nrow, ncol = npar)
  col_index = 1
  
  # add columns for variables with generic cx
  for(com in generic){	
    X[,col_index] = as.vector(t(as.matrix(D[,com])))
    col_index = col_index + 1
  }
  
  # add columns for variables with specific cx
  line_index = 1
  for(spec in specific){
    
    for(v in spec){
      X[seq(from = line_index, by = nalt, length = n),col_index] = D[,v] 
      col_index = col_index + 1
    }
    line_index = line_index + 1
  }
  
  labels = c()
  for(var_group in generic)
    labels = c(labels,paste(var_group,collapse="_"))
  for(var_group in specific)
    labels = c(labels,var_group)
  colnames(X) = labels
  
  return(X)
}

#' function to get utility based model design matrix
#' 
#' \code{getDiscreteMatrix} works like \code{\link{create_X}} but takes specific and generic in a list
#' 
#' @param m a list that contain generic and specific
#' @param D the data set
#' @return the design matrix X
# @export
#getDiscreteMatrix = function(m, D){
#  create_X(m$generic, m$specific, D)
#}

#' calculate the number of variables in a utility specification
#' 
#' \code{getNumVar} calculates the number of columns in design matrix specific in m
#' 
#' @param m a list that contains generic and specific.
#' @param D the data set
#' @return number of parameters
#' @export
getNumVar = function(generic, specific, D){
  ncol(create_X(generic, specific, D))
}

#' formats the names of parameters of utility based model
#' 
#' \code{getNames} computes the names of corresponding parameters in a utility based discrete choice model
#' 
#' @param generic generic parameters
#' @param specific specific parameters
#' @param D the data set
#' @return a vector of the names of each parameter
#' @export
getNames = function(generic, specific, D, ASC = FALSE){
  colnames(create_X(generic, specific, D, ASC))
}

#' adds ASC to specific parameters
augmentSpecific = function(specific, ASC){
	ASCInclude = suppressWarnings(formatASC(ASC, length(specific)))
	for(i in ASCInclude)
		specific[[i]] = c("const",specific[[i]])
	specific
}

formatASC = function(ASC, nAlt){
	ASCInclude = 2:nAlt
	if(1 == nAlt)
		ASCInclude = NULL
	if( is.logical(ASC) & FALSE == ASC)
		ASCInclude = NULL
	if(is.numeric(ASC))
		ASCInclude = ASC
	ASCInclude
}
