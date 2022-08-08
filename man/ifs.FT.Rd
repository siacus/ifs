% file ifs/man/ifs.FT.Rd
% copyright (C) 2001-5 S.M. Iacus    
%
\name{ifs.FT}
\alias{ifs.FT}
\alias{ifs.setup.FT}
\alias{ifs.pf.FT}
\alias{ifs.df.FT}
\alias{IFS.pf.FT}
\alias{IFS.df.FT}
\title{
IFS estimator
}
\description{
  Distribution function estimator based on inverse Fourier transform of ans IFSs.  
}
\usage{
ifs.FT(x, p, s, a, k = 2)
ifs.setup.FT(m, p, s, a, k = 2, cutoff)
ifs.pf.FT(x,b,nterms)
ifs.df.FT(x,b,nterms)
IFS.pf.FT(y, k = 2, n = 512, maps=c("quantile","wl1","wl2")) 
IFS.df.FT(y, k = 2, n = 512, maps=c("quantile","wl1","wl2")) 
}
\arguments{
  \item{x}{where to estimate the function}
  \item{p}{the vector of coefficients \eqn{p_i}}
  \item{s}{the vector of coefficients \eqn{s_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{a}{the vector of coefficients \eqn{a_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{m}{the vector of sample moments}
  \item{k}{number of iterations, default = 2}
  \item{y}{a vector of sample observations}
  \item{n}{the number of points in which to calculate the estimator}
  \item{maps}{type of affine maps}
  \item{b}{the Fourier coefficients}
  \item{nterms}{the number of significant Fourier coefficients after the cutoff}
  \item{cutoff}{cutoff used to determine how many Fourier coefficients are needed}
}
\value{
  The estimated value of the Fourier transform  for \code{ifs.FT}, the estimated value
  of the distribution function for \code{ifs.pf.FT} and the estimated value
  of the density function for \code{ifs.df.FT}. 
  A list of `x' and `y' coordinates plus the Fourier coefficients and the number of
  significant coefficients of the distribution function estimator for \code{IFS.pf.FT} 
  and the density function for \code{IFS.df.FT}.
  The function \code{ifs.setup.FT} return a list of Fourier coefficients and the number
  of significant coefficients.
}
\details{
  This estimator is intended to estimate the continuous distribution
  function, the charateristic function (Fourier transform) and the density function
  of a random variable on [0,1]. 
 }
\note{
 Details of this tecnique can be found in Iacus and La Torre, 2002.
 }
\author{S. M. Iacus}
\references{
  Iacus, S.M, La Torre, D. (2005)
  Approximating distribution functions by iterated function systems, 
  \emph{Journal of Applied Mathematics and Decision Sciences},
  1, 33-46.
}
\seealso{
\code{\link{ecdf}}
}
\examples{
require(ifs)


nobs <- 100
y<-rbeta(nobs,2,4)

# uncomment if you want to test the normal distribution
# y<-sort(rnorm(nobs,3,1))/6


IFS.est <- IFS(y)
xx <- IFS.est$x
tt <- IFS.est$y

ss <- pbeta(xx,2,4)

# uncomment if you want to test the normal distribution   
# ss <- pnorm(6*xx-3)
     
par(mfrow=c(3,1))   
  
plot(ecdf(y),xlim=c(0,1),main="IFS estimator versus EDF")
lines(xx,ss,col="blue")
lines(IFS.est,col="red")
IFS.FT.est <- IFS.pf.FT(y)
xxx <- IFS.FT.est$x
uuu <- IFS.FT.est$y
sss <- pbeta(xxx,2,4)
# uncomment if you want to test the normal distribution   
# sss <- pnorm(6*xxx-3)

lines(IFS.FT.est,col="green")


# calculates MSE


ww <- ecdf(y)(xx)
mean((ww-ss)^2)
mean((tt-ss)^2)
mean((uuu-sss)^2)

plot(xx,(ww-ss)^2,main="MSE",type="l",xlab="x",ylab="MSE(x)")
lines(xx,(tt-ss)^2,col="red")
lines(xxx,(uuu-sss)^2,col="green")

plot(IFS.df.FT(y),type="l",col="green",ylim=c(0,3),main="IFS vs Kernel")
lines(density(y),col="blue")
curve(dbeta(x,2,4),0,1,add=TRUE)
# uncomment if you want to test the normal distribution   
# curve(6*dnorm(x*6-3,0,1),0,1,add=TRUE)
}
\keyword{distribution}

