\name{calRcyRcp}
\alias{calRcyRcp}
\alias{calRcyRcp.default}

\title{
Recycling ratio and recuperation calculation
}

\description{
Calculating recycling ratio, recuperation, and associated standard errors.
}

\usage{
calRcyRcp(Curvedata, Ltx)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): a three-column matrix (i.e., regenerative doses,\cr
 sensitivity-corrected regenerative-dose signals, and associated standard errors)}

  \item{Ltx}{\link{vector}(\bold{required}): a two-element vector consists of sensitivity-corrected\cr
natural-dose signal and its error}
}


\value{
 Return a \link{list} that contains the following elements:
  \item{RecyclingRatio1}{the first recycling ratio and its standard error}
  \item{RecyclingRatio2}{the second recycling ratio and its standard error}
  \item{RecyclingRatio3}{the third recycling ratio and its standard error}
  \item{Recuperation1}{the first recuperation (i.e., ratio of the sensitivity-corrected zero-dose signal to natural-dose signal) 
and its standard error in percent}
  \item{Recuperation2}{the second recuperation (i.e., ratio of the sensitivity-corrected zero-dose signal to maximum regenerative-dose signal) 
and its standard error in percent}
}

\references{
Wintle AG, Murray AS, 2006. A review of quartz optically stimulated luminescence characteristics and
their relevance in single-aliquot regeneration dating protocols. Radiation Measurements, 41(4): 369-391.
}

\seealso{
  \link{calED}; \link{fastED};
  \link{calSARED}; \link{pickSARdata}
}

\note{
If the sensitivity-corrected signals for the frist, second, and third repeated regenerative doses are R1, R2, and R3, respectively, 
then \code{RecyclingRatio1=R2/R1}, \code{RecyclingRatio2=R3/R1}, and \cr
\code{RecyclingRatio3=R3/R2}.
}

\keyword{growth curve}
