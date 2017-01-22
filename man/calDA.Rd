\name{calDA}
\alias{calDA}
\alias{calDA.default}
\title{
Annual dose rate (age) calculation
}
\description{
Calculating the annual dose rate and burial age.
}
\usage{
calDA(dose, minGrainSize, maxGrainSize, Ucontent, Thcontent, Kcontent, 
      Wct, depth, altitude, latitude, longitude, bulkDensity = 2.5, 
      alphaValue = 0.03, nsim = 10000, rdcf = 0.05, rba = 0.05, 
      plot = TRUE)
}

\arguments{
  \item{dose}{\link{vector}(\bold{required}): equivalent dose and associated measurement error (unit, Gy)}

  \item{minGrainSize}{\link{numeric}(\bold{required}): lower limit on grain size (unit, um)}

  \item{maxGrainSize}{\link{numeric}(\bold{required}): upper limit on grain size (unit, um)}

  \item{Ucontent}{\link{vector}(\bold{required}): uranium content and its measurement error (unit, ppm)}

  \item{Thcontent}{\link{vector}(\bold{required}): thorium content and its measurement error (unit, ppm)}

  \item{Kcontent}{\link{vector}(\bold{required}): potassium content and its measurement error (unit, 1 percent)}

  \item{Wct}{\link{vector}(\bold{required}): water content and its measurement error (unit, 1 percent)}

  \item{depth}{\link{numeric}(\bold{required}): sampling depth (unit, m)}

  \item{altitude}{\link{numeric}(\bold{required}): altitude of the sampling site (unit, m above sea level)}

  \item{latitude}{\link{numeric}(\bold{required}): latitude of the sampling site (unit, decimal degree)}

  \item{longitude}{\link{numeric}(\bold{required}): longitude of the sampling site (unit, decimal degree)}

  \item{bulkDensity}{\link{numeric}(with default): average density of bulk sample (unit, g/cm^3)}

  \item{alphaValue}{\link{numeric}(with default): average alpha efficiency}

  \item{nsim}{\link{integer}(with default): number of Monte Carlo iterations}

  \item{rdcf}{\link{numeric}(with default): constant relative standard error for dose-rate conversion factors (unit, 1 percent)}

  \item{rba}{\link{numeric}(with default): constant relative standard error for beta dose absorption fraction (unit, 1 percent)}

  \item{plot}{\link{logical}(with default): draw a plot or not}
}

\details{
Function \link{calDA} is used for calculating the annual dose rate and burial age using concentrations of radioactive nuclides (uranium, thorium, potassium) obtained from Neutron Activation Analysis (NAA), grain size, water content, average sample density, geographical parameters (depth, altitude, latitude, longitude), and an equivalent dose value. The elemental concentrations are converted into annual dose rate according to dose-rate conversion factors (Adamiec and Aitken, 1998). Beta dose absorded fractions are calculated through linear interpolation using published data (Fain et al., 1999). The cosmic ray dose rate is estimated as a function of depth, altitude and geomagnetic latitude (Prescott and Hutton, 1994).\cr\cr
The standard error of an annual dose rate is estimated by a "parametric bootstrap" method. To perform the simulation, dose-rate conversion factors and beta dose absorption factor are assumed to have constant relative standard errors. Arguments such as \code{dose}, \code{Ucontent}, \code{Thcontent}, \code{Kcontent}, \code{wct} should be two-element vectors with the form \code{c(value,sd)}, as their measurement uncertainties are major error sources. Arguments such as \code{depth}, \code{altitude}, \code{latitude}, \code{longitude}, \code{bulkDensity}, \code{alphaValue} can be either a scalar of a two-element vector. This means that uncertainties from these quantities can be either ignored or taken into account during the simulation.
}

\value{
Return a \link{matrix} that contains calculated annual dose rate and age
}

\references{
Adamiec G, Aitken M, 1998. Dose-rate conversion factors: update. Ancient TL, 16(2): 37-49.

Fain J, Soumana S, Montret M, Miallier D, Pilleyre T, Sanzelle S, 1999. Luminescence and ESR dating Beta-dose attenuation for various grain shapes calculated by a Monte-Carlo method. Quaternary Science Reviews, 18(2): 231-234.

Prescott JR, Hutton JT, 1994. Cosmic ray contributions to dose rates for Luminescence and Esr dating: large depths and long-term time variations. Radiation Measurements, 23(2-3): 497-500.
}

\examples{
calDA(dose=c(25.04,0.68), minGrainSize=90,
      maxGrainSize=180, Ucontent=c(2.86,0.19),
      Thcontent=c(8.63,0.34), Kcontent=c(2.00,0.11),
      Wct=c(0.05,0.05), depth=c(1.1,0.05), altitude=c(1170,58.5), 
      latitude=c(37.64,1.88), longitude=c(103.16,5.16),
      bulkDensity=c(2.5,0.1), nsim=10000, rdcf=0.05, rba=0.05) 
}
\keyword{dose rate}
