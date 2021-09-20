#’ autoRegCoef_design_Toeplitz
#'
#' Generate a toeplitz cholesky factor design matrix
#'
#' Generate design matrix for W, the autoregression coefficient model design matrix in Joint Mean Covariance Model
#' The result is a d*(d-1)/2*d design matrix, multiply the parameter coefficient lambda, then can have a wide expand vector,
#' use BayesJMCM2::vec2mat can transfer back to a lower-triangular matrix,
#' @param d The dimension of the matrix
#' @return The design matrix
#'
#'@export
#'
#' @examples
#' BayesJMCM2::vector2ltril(autoRegCoef_design_Toeplitz(5)%*% c(4,3,2,1))
#' \dontrun{vector2ltril(autoRegCoef_design_Toeplitz(5)%*% c(4,3,2,1))}
autoRegCoef_design_Toeplitz=function(d)
{
  # Toeplitz structure indicate the number of parameters need to be estimated is d-1
  # w8! this design is the toeplitz for cholesky factor, not the covariance matrix itself, need more check.
  # Yeah, but this could not be a problem, since old (ti-tj)^q is also have the "Toeplitz" structure in Bayes factor,
  # this idea is also work, but not that interesting

  # This function generate a design matrix W corresponding to Toeplitz structure after vec2Mat
  Des_mat=matrix(rep(0,(d-1)*(d)/2*(d-1)),ncol=d-1)
  count=1
  for (i in 2:d)
  {
    for (j in 1:(i-1))
    {
      Des_mat[count,abs(i-j)]=1
      count=count+1
    }
  }
  return(Des_mat)
}






