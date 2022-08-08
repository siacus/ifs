% file ifs/man/IFSM.Rd
% copyright (C) 2005 S.M. Iacus    
%
\name{IFSM}
\alias{IFSM}
\title{
IFSM operator
}
\description{
  IFSM operator
}
\usage{
IFSM(x, cf, a, s, k = 2)
}
\arguments{
  \item{x}{where to approximate the function}
  \item{cf}{the vector of coefficients \eqn{phi_i}}
  \item{s}{the vector of coefficients \eqn{s_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{a}{the vector of coefficients \eqn{a_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{k}{number of iterations, default = 2}
}
\value{
  The value of the approximate target function.
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
\examples{
require(ifs)

set.seed(123)
n <- 50
dt <- 1/n
t <- (1:n)*dt
Z <- rnorm(n)
B <- sqrt(dt)*cumsum(Z)

ifsm.w.maps() -> maps
a <- maps$a
s <- maps$s

ifsm.setQF(B, s, a) -> QF
ifsm.cf(QF$Q,QF$b,QF$L1,QF$L2,s)-> SOL
psi <- SOL$psi

t1 <- seq(0,1,length=250)
as.numeric(sapply(t1, function(x) IFSM(x,psi,a,s,k=5))) -> B.ifsm
old.mar <- par()$mar
old.mfrow <- par()$mfrow
par(mfrow=c(2,1))
par(mar=c(4,4,1,1))
plot(t1,B.ifsm,type="l",xlab="time",ylab="IFSM")
plot(t,B,col="red",type="l",xlab="time",ylab="Euler scheme")
par(mar=old.mar)
par(mfrow=old.mfrow)
}
\keyword{misc}

