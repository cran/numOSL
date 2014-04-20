\name{dbED}
\alias{dbED}
\alias{dbED.default}
\title{
Summarizing the distribution of a number of equivalent dose values
}
\description{
Calculating statistical parameters (mean, standard deviation, skewness, kurtosis, over-dispersion); performing a homogeneity test; drawing a histogram and a probability density plot.
}
\usage{
dbED(EDdata, aliquot=c("sa","sg"), plot = TRUE, from = NULL, 
     to = NULL, step = 0.1, nbin = 15, samplename = NULL)
}
\arguments{
  \item{EDdata}{\link{data.frame}(\bold{required}): a number of equivalent doses and their standard errors (two columns)}
  \item{aliquot}{\link{character}(with default): single-aliquot ("sa") or single-grain ("sg"), default aliquot="sa"}
  \item{plot}{\link{logical}(with default): draw a histogram and a probability density plot or not}
  \item{from}{\link{numeric}(optional): the desired lower limit on x-axis}
  \item{to}{\link{numeric}(optional): the desired upper limit on x-axis}
  \item{step}{\link{numeric}(with default): a step size used for constructing the probability density plot. Finer step value gives smoother density curve}
  \item{nbin}{\link{numeric}(with default): the desired number of intervals in the histogram}
  \item{samplename}{\link{character}(optional): the name of the sample}
}
\value{
Return a list that contains following elements:
  \item{summaryED}{a basically statistical summary about the given equvalent dose values}
  \item{mlpED}{the mean of the lowest percent [x] of equivalent dose values (Olley et al, 1998)}
  \item{homogeneity}{the results of simple homogeneity test (Galbraith, 2003)}
  \item{skewness}{the skewness of equivalent dose values (Bailey and Arnold, 2006)}
  \item{kurtosis}{the kurtosis of equivalent dose values (Bailey and Arnold, 2006)}
  \item{proposal}{a proposed statistical age model (Bailey and Arnold, 2006; Arnold et al, 2007) that depends on the type of aliquot ("sg" or "sa"), which cannot be relied on absolutely}
}
\references{
Olley, J.M., Caitcheon, G., Murray, A.S., 1998. The distribution of apparent dose as determined by optically stimulated luminescence in small aliquots of fluvial quartz: implications for dating young
sediments. Quaternary Geochronology, 17, pp.1033-1040.

Galbraith, R., 2003. A simple homogeneity test for estimates of dose obtained using OSL. Ancient TL, 21(2), pp.75-77.

Bailey, R.M., Arnold, L.J., 2006. Statistical modelling of single grain quartz De distributions and an assessment of procedures for estimating burial dose. Quaternary Science Reviews, 25 (19-20), pp. 2475-2502.

Arnold, L.J., Bailey, R.M., Tucher, G.E., 2007. Statistical treatment of fluvial dose distributions from southern Colorado arroyo deposits. Quaternary Geochronology, 2, pp.162-167.

Duller, G.A.T., 2008. Single-grain optical dating of Quaternary sediments: why aliquot size matters in luminescence dating. Boreas, 37 (4), pp.589-612.
}
\seealso{
\link{RadialPlotter}
}
\examples{
 data(EDdata)
 dbED(EDdata$al3)
}
\keyword{equivalent dose}
