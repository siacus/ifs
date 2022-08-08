% file ifs/man/ifsm.cf.Rd
% copyright (C) 2005 S.M. Iacus    
%
\name{ifsm.cf}
\alias{ifsm.cf}
\title{Calculates the main parameters of the IFSM operator}
\description{
  Tool function to construct and find the solution of the minimization
  problem involving the quadratic form \eqn{x'Qx + b'x}. Not an optimal one.
  You can provide one better then this.
}
\usage{
ifsm.cf(Q, b, d, l2, s, mu=1e-4)
}
\arguments{
  \item{Q}{the matrix \eqn{Q} of \eqn{x'Qx + b'x}}
  \item{b}{the vector \eqn{b} of \eqn{x'Qx + b'x}}
  \item{d}{the L1 norm of the target function}
  \item{l2}{the L2 norm of the target function}
  \item{s}{the vector \eqn{s} in: \eqn{w_i = s_i *x + a_i}}
  \item{mu}{tolerance}  
}
\value{
	A list
	\item{cf}{the vector of the coefficients to be plugged into the IFSM}
	\item{delta}{the collage distance at the solution}
}
\references{
  Iacus, S.M, La Torre, D. (2005)
  IFSM representation of Brownian motion with applications to simulation, 
  \emph{forthcoming}.
}
\seealso{
\code{\link{IFSM}}
}
\keyword{misc}


