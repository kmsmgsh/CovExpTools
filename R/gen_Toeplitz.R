#' Generate Toeplitz structure matrix given first row of matrix
#'
#' @param toepliz A vector correspond the first row of Toeplitz matrix
#' @return The corresponding toeplitz matrix
#' @examples
#' ToeplitzMat(c(1,2,3))
#' ToeplitzMat(c(1,0.5,0.5))
#' ToeplitzMat(cumprod(c(1,rep(0.7,5))))
#' ToeplitzMat(c(1,1,0,0))
#' @export
ToeplitzMat=function(toeplitz)
{
  p=length(toeplitz)
  A=matrix(rep(0,p^2),ncol=p)
  for (i in 1:p)
  {
    A[i,i:p]=toeplitz[1:(p-i+1)]
  }
  C=t(A)
  C[upper.tri(A)]=A[upper.tri(A)]
  return(C)
}


#' Generate AR1 structure matrix
#' @param rho a numeric value, \eqn{\rho} in AR structure
#' @param p a integer, number of column of matrix
#' @param band=NULL a integer, banded AR(1) covariance matrix, default is not banded.
#' @export
AR1_mat=function(rho,p,band=NULL)
{
  if (is.null(band)){
  toep=cumprod(c(1,rep(rho,p-1)))
  }
  else
  {
    if (band>p)
      error("band out of bound")
    else
    {
      toep=c(cumprod(c(1,rep(rho,band))),rep(0,p-1-band))
    }
  }
  return(ToeplitzMat(toep))
}

