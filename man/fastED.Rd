\name{fastED}
\alias{fastED}
\alias{fastED.default}
\title{
Fast-component equivalent dose calculation
}
\description{
Estimating a fast-component equivalent dose using decay curves obtained from the single aliquot regenerative-dose (SAR) method.
}
\usage{
fastED(Sigdata, Redose, delay.off = c(0,0), ncomp = 2, 
       constant = TRUE, control.args = list(), typ = "cw", 
       model = "gok", origin = FALSE, nsim = 1000, 
       weight.decomp = FALSE, weight.fitGrowth = TRUE,
       trial = TRUE, outpdf = NULL, log = "x", lwd = 2, 
       test.dose = NULL, agID = NULL)
}
\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a series of decay curves stored in a matrix column by column, the first column denotes stimulation time values, see details. 
Data structure of this kind can be obtained using function \link{pickBINdata} by setting argument \code{force.matrix=TRUE}, see examples}

  \item{Redose}{\link{vector}(\bold{required}): regenerative dose values. Example: \code{Redose=c(1,2,3,4,0,1)}}

  \item{delay.off}{\link{vector}(with default): a two-elment vector indicating the "Delay" and "Off"\cr
values of the decay curves, i.e., \code{delay.off[1]=Delay},\code{delay.off[2]=Off}}

  \item{ncomp}{\link{integer}(with default): number of decomposed components}

  \item{constant}{\link{logical}(with default): logical value indicating if a constant background should be subtracted from the decay curve, see function \link{decomp} for details}

  \item{control.args}{\link{list}(with default): arguments used in the differential evolution algorithm, see function \link{decomp} for details}

  \item{typ}{\link{character}(with default): type of an OSL decay curve, only CW-OSL decay curve can be analyzed currently}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see function\cr 
\link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): logical value indicating if the growth curve should be forced to pass the origin}

  \item{nsim}{\link{integer}(with default): number of Monte Carlo iterations used for equivalent dose error estimation}

  \item{weight.decomp}{\link{character}(with default): logical value indicating if the decay curve should be fitted using a weighted procedure, see function \link{decomp} for details}

  \item{weight.fitGrowth}{\link{character}(with default): logical value indicating if the growth curve should be fitted using a weighted procedure, see function \link{fitGrowth} for details}

  \item{trial}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using other models if the given model fails, see function \link{fitGrowth} for details}

  \item{outpdf}{\link{character}(optional): if specified, results of fast-component equivalent dose calculation will be written to a PDF file 
named \code{"outpdf"} and saved to the current work directory}

  \item{log}{\link{character}(with default): a character string which contains "x" if the x axis is to be logarithmic, "y" if the y axis is to be logarithmic and "xy" or "yx" if both axes are to be logarithmic}

  \item{lwd}{\link{numeric}(with default): width of curves (lines)}

  \item{test.dose}{\link{numeric}(optional): test dose of decay curves}

  \item{agID}{\link{vector}(optional): a three-elemenet vector indicating aliquot (grain) ID, i.e.,\cr 
\code{agID[1]=NO}, \code{agID[2]=Position}, \code{agID[3]=Grain}}
}

\details{
Function \link{fastED} is used to estimate a fast-component equivalent dose using data sets obtained from the SAR protocol (Murray and Wintle, 2000). The routine trys to decompose a series of decay curves to a specified number of components, then the numbers of trapped electrons from the fast-component will be used to construct the growth curve to estimate a fast-component equivalent dose. See function \link{decomp}, \link{fitGrowth}, and \link{calED} for more details concerning decay curve decomposition, growth curve fitting, and equivalent dose calculation, respectively.\cr\cr
Argument \code{Sigdata} is a column-matrix made up with stimulation time values and a number of decay curves:
\tabular{lll}{
\bold{Column.no} \tab  \bold{Description}                                     \cr
\code{I}         \tab   Stimulation time values                               \cr
\code{II}        \tab   Natural-dose signal values                            \cr
\code{III}       \tab   Test-dose signal values for the natural-dose          \cr
\code{IV}        \tab   The 1th Regenerative-dose signal values               \cr 
\code{V}         \tab   Test-dose signal values for the 1th regenerative-dose \cr
\code{VI}        \tab   The 2th regenerative-dose signal values               \cr
\code{VII}       \tab   Test-dose signal values for the 2th regenerative-dose \cr
\code{...}       \tab   ...
}
}

\value{
Return an invisible \link{list} containing the following elements:
  \item{decomp.pars}{a \link{list} containing optimized parameters of successfully fitted decay curves}
  \item{Curvedata}{data sets used for building the fast-component growth curve}
  \item{Ltx}{sensitivity-corrected natural-dose fast-component signal and its standard error}
  \item{LMpars}{optimizaed parameters for the fast-component growth curve}
  \item{value}{minimized objective for the fast-component growth curve}
  \item{avg.error}{average fit error for the fast-component growth curve}
  \item{RCS}{reduced chi-square value for the fast-component growth curve}
  \item{FOM}{figure of merit value for the fast-component growth curve in percent}
  \item{calED.method}{method used for fast-component equivalent dose calculation, i.e.,\cr 
\code{"Interpolation"} or \code{"Extrapolation"}}
  \item{mcED}{randomly simulated fast-component equivalent doses}
  \item{ED}{fast-component equivalent dose and its standard error}
  \item{ConfInt}{68 percent and 95 percent confidence interval of fast-component equivalent dose}
  \item{RecyclingRatio1}{the first fast-component recycling ratio and its standard error}
  \item{RecyclingRatio2}{the second fast-component recycling ratio and its standard error}
  \item{RecyclingRatio3}{the third fast-component recycling ratio and its standard error}
  \item{Recuperation1}{the first fast-component recuperation (i.e., ratio of the sensitivity-corrected\cr
zero-dose signal to natural-dose signal) and its standard error in percent}
  \item{Recuperation2}{the second fast-component recuperation (i.e., ratio of the sensitivity-corrected zero-dose 
signal to maximum regenerative-dose signal) and its standard error in percent}

}
\references{
Li SH, Li B, 2006. Dose measurement using the fast component of LM-OSL signals from quartz. 
Radiation Measurements, 41(5): 534-541.

Murray AS, Wintle AG, 2000. Luminescence dating of quartz using improved single-aliquot regenerative-dose protocol. 
Radiation Measurements, 32(1): 57-73.
}

\seealso{
  \link{pickBINdata}; \link{Signaldata};\cr
  \link{fitGrowth}; \link{decomp}; \link{calED}
}

\note{
Argument \code{test.dose} and \code{agID} have nothing to do with fast-component equivalent dose calculation.\cr\cr
The number of trapped electrons that corresponds to the largest decay rate will be regarded as the fast-component signal, which cannot always ensure that a pure fast-component signal be extracted if an ultra-fast decaying component appears.\cr\cr
The authors thank Professor Sheng-Hua Li and Professor Geoff Duller for their helpful discussions concerning fast-component equivalent dose calculation.
}

\examples{
 ### Example 1 (not run):
 # data(Signaldata)
 # fastED(Signaldata$cw,Redose=c(80,160,240,320,0, 80)*0.13,
 #        ncomp=3, constant=FALSE, outpdf="fastED1")

 ### Example 2 (not run):
 # data(BIN)
 # obj_pickBIN <- pickBINdata(BIN, Position=6, Grain=0, 
 #                            LType="OSL", force.matrix=TRUE)
 # fastED(obj_pickBIN$BINdata[[1]], ncomp=2, constant=TRUE,
 #        Redose=c(100,200,300,400,0,100)*0.13, outpdf="fastED2")

}

\keyword{decay curve}
\keyword{growth curve}
\keyword{equivalent dose}
