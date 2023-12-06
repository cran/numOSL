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
      Ucontent, Thcontent, Kcontent, Rbcontent, Wct, depth, longitude, 
      latitude, altitude, alphaValue = 0, inKcontent = 0, inRbcontent = 0, 
      calRbfromK = FALSE, bulkDensity = 2.5, cfType = "Liritzis2013", rdcf = 0, 
      rba = 0, ShallowGamma = TRUE, nsim = 5000, reject = TRUE, plot = TRUE, 
      sampleName = "") 

calDAbatch(inputfile = "inputDRtable", cfType = "Liritzis2013", 
           rdcf = 0, rba = 0, calRbfromK = FALSE,  
           ShallowGamma = TRUE,  nsim = 5000, reject = TRUE, 
           outfile = paste(inputfile,"_Results",sep=""), 
           outpdf = paste(inputfile,"_Results",sep=""), digits = 4)
}

\arguments{
  \item{dose}{\link{vector}(\bold{required}): a two-element vector containing the equivalent dose and associated measurement error (unit, Gy)}

  \item{minGrainSize}{\link{numeric}(\bold{required}): lower limit on grain size (unit, um)}

  \item{maxGrainSize}{\link{numeric}(\bold{required}): upper limit on grain size (unit, um)}

  \item{Ucontent}{\link{vector}(\bold{required}): a two-element vector containing the uranium content and its measurement error (unit, ppm)}

  \item{Thcontent}{\link{vector}(\bold{required}): a two-element vector containing the thorium content and its measurement error (unit, ppm)}

  \item{Kcontent}{\link{vector}(\bold{required}): a two-element vector containing the potassium content and its measurement error (unit, percent)}

  \item{Rbcontent}{\link{numeric}(\bold{required}): the rubidium content (unit, ppm). The measurement error is optional}  

  \item{Wct}{\link{vector}(\bold{required}): a two-element vector containing the water content and its measurement error (unit, percent)}

  \item{depth}{\link{numeric}(\bold{required}): sampling depth (unit, m). The associated error is optional}

  \item{longitude}{\link{numeric}(\bold{required}): longitude of the sampling site (unit, decimal degree). The associated error is optional}

  \item{latitude}{\link{numeric}(\bold{required}): latitude of the sampling site (unit, decimal degree). The associated error is optional}

  \item{altitude}{\link{numeric}(\bold{required}): altitude of the sampling site (unit, m above sea level). The associated error is optional}

  \item{alphaValue}{\link{numeric}(with default): average alpha efficiency. The associated error is optional, for example, \code{alphaValue=0.038} or \code{alphaValue=c(0.038,0.002)}}

  \item{inKcontent}{\link{numeric}(with default): internal potassium content (unit, percent). The associated error is optional, for example, \code{inKcontent=12.5}, or \code{inKcontent=c(12.5,0.5)}}

  \item{inRbcontent}{\link{numeric}(with default): internal rubidium content (unit, ppm). The associated error is optional, for example, \code{inRbcontent=400}, or \code{inRbcontent=c(400,100)}}

  \item{calRbfromK}{\link{logical}(with default): whether calculate the rubidium and internal rubidium contents using the provided potassium and internal potassium contents respectively.  
If \code{calRbfromK=TRUE}, the provided rubidium and/or internal rubidium contents will not be used for dose-rate calculation}

  \item{bulkDensity}{\link{numeric}(with default): average density of bulk sample (unit, g/cm^3).\cr The associated error is optional, for example, \code{bulkDensity=2.5},\cr 
or \code{bulkDensity=c(2.5,0.2)}}

  \item{cfType}{\link{character}(with default): type of the conversion factor, one from\cr \code{"AdamiecAitken1998"}, \code{"Guerin2011"}, and \code{"Liritzis2013"}}

  \item{rdcf}{\link{numeric}(with default): constant relative standard error (RSD) for dose-rate conversion factors (unit, percent). If \code{rdcf=0}, the dose-rate conversion factors will be assummed to be perfectly determined, otherwise their errors calculated using the constant RSD will be accounted for during the Monte Carlo simulation. Since the conversion factors of \code{"Liritzis2013"} contain measured standard errors, when \code{cfType="Liritzis2013"}, a positive \code{rdcf} value indicates that the measured errors will be accounted for during simulation}

  \item{rba}{\link{numeric}(with default): constant RSD for alpha and beta dose absorption fractions (unit, percent). If \code{rba=0}, the determined alpha and beta dose attenuation factors will be assummed to be free from errors}

  \item{ShallowGamma}{\link{logical}(with default): consider the scaling of gamma dose rate for samples collected at shallow burial depths or not}

  \item{nsim}{\link{integer}(with default): number of Monte Carlo simulations}

  \item{reject}{\link{logical}(with default): whether randomly generated negative values of variables (including equivalent dose, uranium, thorium, potassium, and water contents, alpha efficiency, and bulk-sample density, etc) and meaningless values (longitude beyonds \code{[-180,180]}, or latitude beyonds \code{[-90,90]}) will be reject during the Monte Carlo simulation}

  \item{plot}{\link{logical}(with default): draw a plot showing the simulated distributions of dose rates and ages or not}

  \item{sampleName}{\link{character}(with default): name of the sample}

  \item{inputfile}{\link{character}(with default): name of the input file containing measured dataset of individual samples used for dose rate and age calculations}

  \item{outfile}{\link{character}(with default): name of the files containing the results of calculated dose rates and ages for a number of samples. The files will be written using CSV/HTML format and saved to the current work directory}

  \item{outpdf}{\link{character}(with default): name of a PDF file showing the distributions of dose rates and ages simulated using a Monte Carlo method for a number of samples. The file will be saved to the current work directory}

  \item{digits}{\link{integer}(with default):  the number of decimal places or significant digits to be used for values of the output CSV/HTML files}

}

\details{
Function \link{calDA} is used for calculating the annual dose rate and burial age using the concentrations of radioactive nuclides (uranium, thorium, potassium) obtained from Neutron Activation Analysis (NAA) or inductively coupled plasma mass spectrometry (ICP-MS), grain size, water content, average sample density, geographical parameters (depth, longitude, latitude, altitude), and the equivalent dose. The elemental concentrations of uranium, thorium, and potassium are converted into alpha, bata, and gamma dose rates according to dose-rate conversion factors (Adamiec and Aitken, 1998; Guerin et al., 2011; Liritzis et al., 2013). Alpha and beta dose absorded fractions are determined using published data (Mejdahl, 1979; Brennan et al., 1991). The contribution of the internal beta dose rate can be assessed if the internal potassium content is provided. The gamma dose rate for samples collected at shallow depths are corrected using the scaling factors of Aitken (1985). The cosmic ray dose rate is estimated as a function of depth, altitude and geomagnetic latitude (Prescott and Hutton, 1994) and the contribution to cosmic dose rate from a soft component is considered at shallow depths (Barbouti and Rastin, 1983).\cr\cr
The hydrofluoric acid-etched quartz and K-feldspar samples have an alpha efficiency of zero, while the reported alpha values of un-etched coarse-grained quartz and K-feldspar are \eqn{0.1\pm0.02} (Olley et al., 1998) and \eqn{0.15\pm0.05} (Balescu and Lamothe, 1994), respectively. Three reported alpha values for fine-grained quartz are \eqn{0.038\pm0.002} (Rees-Jones, 1995), \eqn{0.04\pm0.01} (Rees-Jones and Tite, 1997), and \eqn{0.035\pm0.003} (Lai et al., 2008). Two reported alpha values for fine-grained polymineral are \eqn{0.086\pm0.004} (Rees-Jones, 1995) and \eqn{0.09\pm0.02} (Kreutzer et al., 2014). Huntley and Hancock (2001) assumed an internal rubidium content of \eqn{400\pm100} ppm. Three reported internal potassium contents are \eqn{12.5\pm0.5\%} (Huntley and Baril, 1997), \eqn{13\pm1\%} (Zhao and Li, 2005), and \eqn{10\pm2\%} (Smedley et al., 2012).\cr\cr 
The standard error of the age and dose rate is estimated by a "parametric bootstrap" method (Peng et al., 2016). Constant relative standard errors for dose-rate conversion factors, alpha and beta dose absorption factors can be assummed during the simulation. Arguments such as \code{dose}, \code{Ucontent}, \code{Thcontent}, \code{Kcontent}, \code{Wct} should be two-element vectors representing paired values of \eqn{\mu\pm\sigma}, where \eqn{\mu} and \eqn{\sigma} denote the measured value and associated standard error, respectively. Arguments such as \code{depth}, \code{longitude}, \code{latitude}, \code{altitude}, \code{alphaValue}, \code{Rbcontent}, \code{inKcontent}, \code{inRbcontent}, and \code{bulkDensity}, can be either a scalar or a two-element vector. This means that uncertainties from these quantities can be either ignored or taken into account during the simulation.\cr\cr
Function \link{calDAbatch} is a wrapper of the function \link{calDA} and is used to calculate the dose rates and burial ages for a number of samples in a batch mode. The function requires as input a CSV file containing dose-rate datasets of different samples that are saved row by row. A template of the input CSV file with the name \code{"myDRdata"} can be generated using the command \code{calDAbatch("myDRdata")}\cr(see examples).
}

\value{
Function \link{calDA} returns a \link{matrix} that contains calculated alpha, beta, gamma, cosmic, and total dose rate and age and associated standard errors, lower and upper bounds of 95 percent confidence intervals for the sample under considered.\cr
Function \link{calDAbatch} returns an invisible \link{list} that contains calculated dose-rate results for each of the provided samples. 
}

\references{
Adamiec G, Aitken M, 1998. Dose-rate conversion factors: update. Ancient TL, 16(2): 37-49.

Aitken MJ, 1985. Thermoluminescence Dating. Academic Press, London.

Balescu S, Lamothe M, 1994. Comparison of TL and IRSL age estimates of feldspar coarse grains from waterlain sediments. Quaternary Science Reviews, 13: 437-444.

Barbouti AI, Rastin BC, 1983. A study of the absolute intensity of muons at sea level and under various thicknesses of absorber. Journal of Physics G Nuclear Physics, 9: 1577e1595.

Brennan BJ, Lyons RG, Phillips SW, 1991. Attenuation of alpha particle track dose for spherical grains. International Journal of Radiation Application and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 18: 249-253.

Guerin G, Mercier N, Adamiec G, 2011. Dose-rate conversion factors: update. Ancient TL, 29: 5-8.

Huntley DJ, Baril M, 1997. The K content of the K-feldspars being measured in optical dating or in thermoluminescence dating. Ancient TL, 15: 11-13.

Huntley DJ, Hancock R, 2001. The Rb contents of the K-feldspar grains being measured in optical dating. Ancient TL, 19: 43-46.

Kreutzer S, Schmidt C, DeWitt R, Fuchs M, 2014. The a-value of polymineral fine grain samples measured with the post-IR IRSL protocol. Radiation Measurements, 69: 18-29.

Lai ZP, Zoller L, Fuchs M, Bruckner H, 2008. Alpha efficiency determination for OSL of quartz extracted from Chinese loess. Radiation Measurements, 43: 767-770. 

Liritzis I, Stamoulis K, Papachristodoulou C, Ioannides K, 2013. A re-evaluation of radiation dose-rate conversion factors. Mediterranean Archaeology and Archaeometry, 13: 1-15.

Mejdahl V, 1979. Thermoluminescence dating: beta-dose attenuation in quartz grains. Archaeometry, 21: 61-72.

Olley J, Caitcheon G, Murray A, 1998. The distribution of apparent dose as determined by Optically Stimulated Luminescence in small aliquots of fluvial quartzImplications for dating young sediments. Quaternary Science Reviews, 17: 1033-1040. 

Prescott, JR, Hutton JT, 1994. Cosmic ray contributions to dose rates for Luminescence and Esr dating: large depths and long-term time variations. Radiation Measurements, 23(2-3): 497-500.

Peng J, Dong ZB, Zhang ZC, 2016. Determining the error of dose rate estimates on luminescence dating using Monte Carlo approach. Nuclear Techniques, 38(7): 070201. (In Chinese).

Rees-Jones J, 1995. Optical dating of young sediments using fine-grained quartz. Ancient TL, 13: 9-14.

Rees-Jones J, Tite MS, 1997. Optical dating results for British archaeological sediments. Archaeometry, 39: 177-187.

Smedley RK, Duller GAT, Pearce NJG, Roberts HM, 2012. Determining the K-content of single-grains of feldspar for luminescence dating. Radiation Measurements, 47: 790-796.

Zhao H, Li SH, 2005. Internal dose rate to K-feldspar grains from radioactive elements other than potassium. Radiation Measurements, 40: 84-93.

\bold{Further reading}

Durcan JA, King GE, Duller GAT, 2015. DRAC: Dose Rate and Age Calculator for trapped charge dating. Quaternary Geochronology, 28: 54-61.

Grun R, 2009. The "AGE" program for the calculation of luminescence age estimates. Ancient TL, 27: 45-46.

}

\author{
Orignal R code written by Jun Peng, improved version of code written by Chunxin Wang
}

\examples{
calDA(dose=c(25.04,0.68), minGrainSize=90,
      maxGrainSize=180, Ucontent=c(2.86,0.19),
      Thcontent=c(8.63,0.34), Kcontent=c(2.00,0.11), Rbcontent=0,
      Wct=c(0.05,0.05), depth=c(1.1,0.05), longitude=c(103.16,0.1),
      latitude=c(37.64,0.1), altitude=c(1170,58.5), bulkDensity=c(2.5,0.1), 
      rdcf=0.03, rba=0.03) 

# Not run.
# Generate a template of input CSV file "myDRdata" using the following command.
# calDAbatch(inputfile="myDRdata")

# Put your dose rate dataset into the above generated template file "myDRdata.csv", then run 
# the following command to calculate dose rates and ages for your samples in a batch mode.
# d <- calDAbatch(inputfile="myDRdata", nsim=3000)
# idx <- sapply(d, is.matrix)
# h <- t(sapply(d[idx],function(x) c(t(x[6:7,1:2]))))
# colnames(h) <- c("DR","Se.DR","Age","Se.Age")
# print(h)
}
\keyword{Annual dose rate}
\keyword{Alpha dose rate}
\keyword{Beta dose rate}
\keyword{Gamma dose rate}
\keyword{Cosmic dose rate}
