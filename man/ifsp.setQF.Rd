% file ifs/man/ifsp.setQF.Rd
% copyright (C) 2001-5 S.M. Iacus    
%
\name{ifsp.setQF}
\alias{ifsp.setQF}
\title{Sets up the quadratic form for the IFSP}
\description{
  Tool function to construct the quadratic form \eqn{x'Qx + b'x} to be minimized
  to construct the IFSP operator.
}
\usage{
ifsp.setQF(m, s, a, n = 10)
}
\arguments{
  \item{m}{the vector of the sample or true moments of the target function}
  \item{s}{the vector of coefficients \eqn{s_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{a}{the vector of coefficients \eqn{a_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{n}{number of parameter to use in the IFSP operator, default = 10}
}
\value{
	\item{Q}{the matrix of the quadratic form}
	\item{b}{the matrix of the quadratic form}	
}
\details{
  This operator is intended to approximate a continuous distribution
  function of a random variable on [0,1]. If moments are estimated on
  a random sample, then the IFSP operator is an estimator of the
  distribution function of the data.
 }
\author{S. M. Iacus}
\references{
  Iacus, S.M, La Torre, D. (2005)
  Approximating distribution functions by iterated function systems, 
  \emph{Journal of Applied Mathematics and Decision Sciences},
  1, 33-46.
}
\seealso{
\code{\link{ifs}}
}
\keyword{misc}

