% file ifs/man/ifs.Rd
% copyright (C) 2001-5 S.M. Iacus    
%
\name{ifs}
\alias{ifs}
\alias{ifs.flex}
\alias{IFS}
\title{
IFS estimator
}
\description{
  Distribution function estimator based on sample quantiles.  
}
\usage{
ifs(x, p, s, a, k = 5)
ifs.flex(x, p, s, a, k = 5, f = NULL)
IFS(y, k = 5, q = 0.5, f = NULL, n = 512, maps = c("quantile", 
    "wl1", "wl2"))
}
\arguments{
  \item{x}{where to estimate the distribution function}
  \item{p}{the vector of coefficients \eqn{p_i}}
  \item{s}{the vector of coefficients \eqn{s_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{a}{the vector of coefficients \eqn{a_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{k}{number of iterations, default = 5}
  \item{y}{a vector of sample observations}
  \item{q}{the proportion of quantiles to use in the construction of the 
  estimator, default = 0.5. The number of quantiles is the 
  \code{q * length(y)}.} 
  \item{f}{the starting point in the space of distribution functions}
  \item{n}{the number of points in which to calculate the IFS}
  \item{maps}{type of affine maps}
}
\value{
  The estimated value of the distribution function for \code{ifs} and \code{ifs.flex} or
  a list of `x' and `y' coordinates of the IFS(x) graph for \code{IFS}.
}
\details{
  This estimator is intended to estimate the continuous distribution
  function of a random variable on [0,1]. The estimator is a continuous
  function not everywhere differentiable.
 }
\note{
 It is asymptotically as good as the empirical distribution function
  (see Iacus and La Torre, 2001).  
 This function is called by \link{IFS}. If you need to call the function
 several times, you should better use \code{ifs} providing the
 points and coefficients once instead of \code{IFS}.
 Empirical evidence shows that the IFS-estimator is better than the edf (even 
 for very small samples) in the sup-norm metric. It is also better in the MSE
 sense outside of the distribution's tails if the sample quantiles are used
 as points.
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


y<-rbeta(50,.5,.1)

# uncomment if you want to test the normal distribution
# y<-sort(rnorm(50,3,1))/6


IFS.est <- IFS(y)
xx <- IFS.est$x
tt <- IFS.est$y

ss <- pbeta(xx,.5,.1)

# uncomment if you want to test the normal distribution   
# ss <- pnorm(6*xx-3)
     
par(mfrow=c(2,1))   
  
plot(ecdf(y),xlim=c(0,1),main="IFS estimator versus EDF")
lines(xx,ss,col="blue")
lines(xx,tt,col="red")


# calculates MSE


ww <- ecdf(y)(xx)
mean((ww-ss)^2)
mean((tt-ss)^2)

plot(xx,(ww-ss)^2,main="MSE",type="l",xlab="x",ylab="MSE(x)")
lines(xx,(tt-ss)^2,col="red")
}
\keyword{nonparametric}

