\name{calDA}
\alias{calDA}
\alias{calDA.default}
\alias{calDAbatch}
\alias{calDAbatch.default}
\title{
Dose rate and age calculation
}
\description{
Calculating the total dose rate and burial age and assessing associated standard errors using a Monte Carlo method.
}
\usage{
calDA(dose, minGrainSize, maxGrainSize,
      Ucontent, Thcontent, Kcontent, Wct, depth, longitude, 
      latitude, altitude, alphaValue = 0, inKcontent = 0, bulkDensity = 2.5, 
      cfType = "Liritzis2013", rdcf = 0, rba = 0, ShallowGamma = TRUE, 
      nsim = 5000, rejectNeg = TRUE, plot = TRUE, sampleName = "") 

calDAbatch(inputfile = "inputDRtable", cfType = "Liritzis2013", 
           rdcf = 0, rba = 0, ShallowGamma = TRUE, nsim = 5000, 
           rejectNeg = TRUE, outfile = paste(inputfile,"_Results",sep=""), 
           outpdf = paste(inputfile,"_Results",sep=""))
}

\arguments{
  \item{dose}{\link{vector}(\bold{required}): a two-element vector containing the equivalent dose and associated measurement error (unit, Gy)}

  \item{minGrainSize}{\link{numeric}(\bold{required}): lower limit on grain size (unit, um)}

  \item{maxGrainSize}{\link{numeric}(\bold{required}): upper limit on grain size (unit, um)}

  \item{Ucontent}{\link{vector}(\bold{required}): a two-element vector containing the uranium content and its measurement error (unit, ppm)}

  \item{Thcontent}{\link{vector}(\bold{required}): a two-element vector containing the thorium content and its measurement error (unit, ppm)}

  \item{Kcontent}{\link{vector}(\bold{required}): a two-element vector containing the potassium content and its measurement error (unit, percent)}

  \item{Wct}{\link{vector}(\bold{required}): a two-element vector containing the water content and its measurement error (unit, percent)}

  \item{depth}{\link{numeric}(\bold{required}): sampling depth (unit, m). The associated error is optional}

  \item{longitude}{\link{numeric}(\bold{required}): longitude of the sampling site (unit, decimal degree). The associated error is optional}

  \item{latitude}{\link{numeric}(\bold{required}): latitude of the sampling site (unit, decimal degree). The associated error is optional}

  \item{altitude}{\link{numeric}(\bold{required}): altitude of the sampling site (unit, m above sea level). The associated error is optional}

  \item{alphaValue}{\link{numeric}(with default): average alpha efficiency. The associated error is optional, for example, \code{alphaValue=0.038} or \code{alphaValue=c(0.038,0.002)}}

  \item{inKcontent}{\link{numeric}(\bold{required}): internal potassium content (unit, percent). The associated error is optional, for example, \code{inKcontent=12.5}, or \code{inKcontent=c(12.5,0.5)}}

  \item{bulkDensity}{\link{numeric}(with default): average density of bulk sample (unit, g/cm^3).\cr The associated error is optional, for example, \code{bulkDensity=2.5},\cr 
or \code{bulkDensity=c(2.5,0.2)}}

  \item{cfType}{\link{character}(with default): type of the conversion factor, one from\cr \code{"AdamiecAitken1998"}, \code{"Guerin2011"}, and \code{"Liritzis2013"}}

  \item{rdcf}{\link{numeric}(with default): constant relative standard error (RSD) for dose-rate conversion factors (unit, percent). If \code{rdcf=0}, the dose-rate conversion factors will be assummed to be perfectly determined, otherwise their errors calculated using the constant RSD will be accounted for during the Monte Carlo simulation. Since the conversion factors of \code{"Liritzis2013"} contain measured standard errors, when \code{cfType="Liritzis2013"}, a positive \code{rdcf} value indicates that the measured errors will be accounted for during simulation}

  \item{rba}{\link{numeric}(with default): constant RSD for alpha and beta dose absorption fractions (unit, percent). If \code{rba=0}, the determined alpha and beta dose attenuation factors will be assummed to be free from errors}

  \item{ShallowGamma}{\link{logical}(with default): consider the scaling of gamma dose rate for samples collected at shallow burial depths or not}

  \item{nsim}{\link{integer}(with default): number of Monte Carlo simulations}

  \item{rejectNeg}{\link{logical}(with default): whether randomly generated negative values of variables (including uranium, thorium, potassium, and water contents, alpha efficiency, and bulk-sample density, etc) will be reject during the Monte Carlo simulation}

  \item{plot}{\link{logical}(with default): draw a plot showing the simulated distributions of dose rates and ages or not}

  \item{sampleName}{\link{character}(with default): name of the sample}

  \item{inputfile}{\link{character}(with default): name of the input file containing measured dataset of individual samples used for dose rate and age calculations}

  \item{outfile}{\link{character}(with default): name of the file containing the results of calculated dose rates and ages for a number of samples. The file will be written to a CSV file and saved to the current work directory}

  \item{outpdf}{\link{character}(with default): name of a PDF file showing the distributions of dose rates and ages simulated using a Monte Carlo method for a number of samples. The file will be saved to the current work directory}

}

\details{
Function \link{calDA} is used for calculating the annual dose rate and burial age using the concentrations of radioactive nuclides (uranium, thorium, potassium) obtained from Neutron Activation Analysis (NAA) or inductively coupled plasma mass spectrometry (ICP-MS), grain size, water content, average sample density, geographical parameters (depth, longitude, latitude, altitude), and the equivalent dose. The elemental concentrations of uranium, thorium, and potassium are converted into alpha, bata, and gamma dose rates according to dose-rate conversion factors (Adamiec and Aitken, 1998; Guerin et al., 2011; Liritzis et al., 2013). Alpha and beta dose absorded fractions are determined using published data (Mejdahl, 1979; Brennan et al., 1991). The contribution of the internal beta dose rate can be assessed if the internal potassium content is provided. The gamma dose rate for samples collected at shallow depths are corrected using the scaling factors of Aitken (1985). The cosmic ray dose rate is estimated as a function of depth, altitude and geomagnetic latitude (Prescott and Hutton, 1994) and the contribution to cosmic dose rate from a soft component is considered at shallow depths (Barbouti and Rastin, 1983).\cr\cr
The standard error of the total dose rate is estimated by a "parametric bootstrap" method. Constant relative standard errors for dose-rate conversion factors, alpha and beta dose absorption factors can be assummed during the simulation. Arguments such as \code{dose}, \code{Ucontent}, \code{Thcontent}, \code{Kcontent}, \code{Wct} should be two-element vectors with the form \code{c(u,sigma)}, where \code{u} and \code{sigma} denote the measured value and associated standard error, respectively. Arguments such as \code{depth}, \code{longitude}, \code{latitude}, \code{altitude}, \code{alphaValue}, \code{inKcontent}, and \code{bulkDensity}, can be either a scalar or a two-element vector. This means that uncertainties from these quantities can be either ignored or taken into account during the simulation.\cr\cr
Function \link{calDAbatch} is a wrapper of the function \link{calDA} and is used to calculate the dose rates and burial ages for a number of samples in a batch mode. The function requires as input a CSV file containing dose-rate datasets of different samples that are saved row by row. A template of the input CSV file with the name \code{"myDRdata"} can be generated using the command \code{calDAbatch("myDRdata")}\cr(see examples).
}

\value{
Return a \link{matrix} that contains calculated alpha, beta, gamma, cosmic, and total dose rate and age and associated standard errors, lower and upper bounds of 95 percent confidence intervals.
}

\references{
Adamiec G, Aitken M, 1998. Dose-rate conversion factors: update. Ancient TL, 16(2): 37-49.

Aitken, M.J., 1985. Thermoluminescence Dating. Academic Press, London.

Barbouti, A.I., Rastin, B.C., 1983. A study of the absolute intensity of muons at sea level and under various thicknesses of absorber. Journal of Physics G Nuclear Physics, 9: 1577e1595.

Brennan, B.J., Lyons, R.G., Phillips, S.W., 1991. Attenuation of alpha particle track dose for spherical grains. International Journal of Radiation Application and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 18: 249-253.

Guerin, G., Mercier, N., Adamiec, G., 2011. Dose-rate conversion factors: update. Ancient TL, 29: 5-8.

Liritzis, I., Stamoulis, K., Papachristodoulou, C., Ioannides, K., 2013. A re-evaluation of radiation dose-rate conversion factors. Mediterranean Archaeology and Archaeometry, 13: 1-15.

Mejdahl, V., 1979. Thermoluminescence dating: beta-dose attenuation in quartz grains. Archaeometry, 21: 61-72.

Prescott JR, Hutton JT, 1994. Cosmic ray contributions to dose rates for Luminescence and Esr dating: large depths and long-term time variations. Radiation Measurements, 23(2-3): 497-500.

\bold{Further reading}

Durcan, J.A., King, G.E., Duller, G.A.T., 2015. DRAC: Dose Rate and Age Calculator for trapped charge dating. Quaternary Geochronology, 28: 54-61.

Grun, R., 2009. The "AGE" program for the calculation of luminescence age estimates. Ancient TL, 27: 45-46.
}

\author{
Orignal R code written by Jun Peng, improved version of code written by Chunxin Wang
}

\examples{
calDA(dose=c(25.04,0.68), minGrainSize=90,
      maxGrainSize=180, Ucontent=c(2.86,0.19),
      Thcontent=c(8.63,0.34), Kcontent=c(2.00,0.11),
      Wct=c(0.05,0.05), depth=c(1.1,0.05), longitude=c(103.16,0.1),
      latitude=c(37.64,0.1), altitude=c(1170,58.5), bulkDensity=c(2.5,0.1), 
      rdcf=0.03, rba=0.03) 

# Not run.
# Generate a template of input CSV file "myDRdata" using the following command.
# calDAbatch(inputfile="myDRdata")

# Put your dose rate dataset into the above generated template file "myDRdata.csv", then run 
# the following command to calculate dose rates and ages for your samples in a batch mode.
# d <- calDAbatch(inputfile="myDRdata", nsim=6000)
# idx <- sapply(d, is.matrix)
# h <- t(sapply(d[idx],function(x) c(t(x[6:7,1:2]))))
# colnames(h) <- c("DR","Se.DR","Age","Se.Age")
# h
}
\keyword{Annual dose rate}
\keyword{Alpha dose rate}
\keyword{Beta dose rate}
\keyword{Gamma dose rate}
\keyword{Cosmic dose rate}
