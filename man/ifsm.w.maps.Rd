% file ifs/man/ifsm.w.maps.Rd
% copyright (C) 2001-5 S.M. Iacus    
%
\name{ifsm.w.maps}
\alias{ifsm.w.maps}
\title{Set up the parameters for the maps of the IFSM operator}
\description{
 This is called before calling \code{ifsm.setQF} to prepare the
 parameters to be passed in \code{ifsm.setQF}.
}
\usage{
ifsm.w.maps(M=8)
}
\arguments{
  \item{M}{M is such that \code{sum(2^(1:M))} maps are created}
}
\value{
	A list of 
	\item{a}{the vector of the coefficents `a' in the maps}	
	\item{s}{the vector of the coefficents `s' in the maps}	
}
\author{S. M. Iacus}
\seealso{
\code{\link{IFSM}}
}
\keyword{misc}


