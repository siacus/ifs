% file ifs/man/ifsp.cf.Rd
% copyright (C) 2001-5 S.M. Iacus    
%
\name{ifsp.cf}
\alias{ifsp.cf}
\title{Calculates the main parameters of the IFS estimators}
\description{
  Tool function to construct and find the solution of the minimization
  problem involving the quadratic form \eqn{x'Qx + b'x}. Not an optimal one.
  You can provide one better then this.
}
\usage{
ifsp.cf(Q,b)
}
\arguments{
  \item{Q}{the matrix \eqn{Q} of \eqn{x'Qx + b'x}}
  \item{b}{the vector \eqn{b} of \eqn{x'Qx + b'x}}
}
\value{
	\item{p}{the vector of the coefficients to be plugged into the IFS}
}
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


