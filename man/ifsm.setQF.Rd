% file ifs/man/ifsm.setQF.Rd
% copyright (C) 2005 S.M. Iacus    
%
\name{ifsm.setQF}
\alias{ifsm.setQF}
\title{Sets up the quadratic form for the IFSM}
\description{
  Tool function to construct the quadratic form \eqn{x'Qx + b'x + l2} to be 
  minimized under some constraint depending on l1. This is used to construct
  the IFSM operator.
}
\usage{
ifsm.setQF(u, s, a)
}
\arguments{
  \item{u}{the vector of values of the target function u}
  \item{s}{the vector of coefficients \eqn{s_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{a}{the vector of coefficients \eqn{a_i} in: \eqn{w_i = s_i *x + a_i}}
}
\value{
	List of elements
	\item{Q}{the matrix of the quadratic form}
	\item{b}{the matrix of the quadratic form}
	\item{L1}{the L1 norm of the target function}
	\item{L2}{the L2 norm of the target function}
	\item{M1}{the integral of the target function}		
}
\details{
  This operator is intended to approximate a function on L2[0,1].
  If `u' is simulated, then the IFSM can be used to simulate a IFSM
  version of `u'.
   }
\author{S. M. Iacus}
\references{
  Iacus, S.M, La Torre, D. (2005)
  IFSM representation of Brownian motion with applications to simulation, 
  \emph{forthcoming}.
}
\seealso{
\code{\link{IFSM}}
}
\keyword{misc}

