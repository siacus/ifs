% file ifs/man/ifsp.w.maps.Rd
% copyright (C) 2001-5 S.M. Iacus    
%
\name{ifsp.w.maps}
\alias{ifsp.w.maps}
\title{Set up the parameters for the maps of the IFSP operator}
\description{
 This is called before calling \code{ifsp.setQF} to prepare the
 parameters to be passed in \code{ifsp.setQF}.
}
\usage{
ifsp.w.maps(y, maps = c("quantile","wl1","wl2"), qtl)
}
\arguments{
  \item{y}{the vector of the sample observations}
  \item{maps}{type of maps: quantile, wl1 or wl2}
  \item{qtl}{instead of passing the data y you can pass a vector
   of quantiles}
}
\value{
	\item{m}{the vector of the empirical moments}
	\item{a}{the vector of the coefficents `a' in the maps}	
	\item{s}{the vector of the coefficents `s' in the maps}	
	\item{n}{the number of maps}	
}
\author{S. M. Iacus}
\seealso{
\code{\link{ifs}}
}
\keyword{misc}


