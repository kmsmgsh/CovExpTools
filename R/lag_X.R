#' Generate lag covariate from design matrix
#' The resulting design matrix W correspond to column-wise expansion of Phi
#' @param X, design matrix for mean vector
#' @return W, design matrix for autocoefficient coefficient
#' @export
#’
lag_X=function(X,ncolW=NULL)
{
  if (is.null(dim(X)))
    error("X is not a matrix")
  # The for loop should be column wise, i.e. column first, inner for loop is from 1:(i-1), i is current row number,
  # then the outter for loop is row wise, from 2:m_i, the number of observations.
  #ncol=dim(X)[2]
  if (is.null(ncolW))
  {
    ncolW=ncol(X)
  }
  if (ncol(X)<ncolW)
    error("ncol of W is greater than X")
  p=dim(X)[1]
  W=matrix(rep(0,p*(p-1)*ncolW/2),ncol=ncolW)
  count=1
  for (i in 2:nrow(X))
  {
    for(j in 1:(i-1))
    {
      W[count,1:ncolW]=X[i,1:ncolW]-X[j,1:ncolW]
      count=count+1
    }
  }
  return (W)
}
