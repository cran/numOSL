\name{calED}
\alias{calED}
\alias{calED.default}

\title{
Equivalent dose calculation and error assessment
}

\description{
Calculating an equivalent dose and assessing its standard error.
}

\usage{
calED(Curvedata, Ltx, model = "gok", origin = FALSE, 
      errMethod = "sp", nsim = 500, weight = TRUE,  
      trial = FALSE, plot = TRUE, nofit.rgd = NULL, 
      agID = NULL, Tn = NULL, Tn3BG = NULL, 
      TnBG.ratio = NULL, rseTn = NULL, FR = NULL, 
      LnTn.curve = NULL, TxTn = NULL)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): a three-column matrix (i.e., regenerative doses,\cr
 sensitivity-corrected regenerative-dose signals, and associated standard errors)}

  \item{Ltx}{\link{vector}(\bold{required}): a two-element vector consists of sensitivity-corrected\cr
natural-dose signal and its error}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): logical value indicating if the growth curve should be forced to pass the origin}

  \item{errMethod}{\link{character}(with default): method used for equivalent dose error assessment.\cr
\code{"sp"} and \code{"mc"} denote error estimation using the Simple Transformation and Monte Carlo methods, respectively}

  \item{nsim}{\link{integer}(with default): desired number of randomly simulated equivalent dose obtained by Monte Carlo simulation}

  \item{weight}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using a weighted procedure, see function \link{fitGrowth} for details}

  \item{trial}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using other models if the given model fails, see function \link{fitGrowth} for details}

  \item{plot}{\link{logical}(with default): logical value indicating if the results should be plotted}

  \item{nofit.rgd}{\link{integer}(optional): regenerative doses that will not be used during the fitting. 
For example, if \code{nofit.rgd=1} then the first regenerative dose will not be used during growth curve fitting}
 
  \item{agID}{\link{vector}(optional): a three-elemenet vector indicating aliquot (grain) ID, i.e.,\cr
\code{agID[1]=NO}, \code{agID[2]=Position}, \code{agID[3]=Grain}}

  \item{Tn}{\link{vector}(optional): a two-element vector containing value and\cr
standard error of Tn}  

  \item{Tn3BG}{\link{numeric}(optional): 0-1 value indicating if Tn is more than 3 sigma above BG,\cr
1 indicates Tn>3_sigma_BG, 0 indicates Tn<=3_sigma_BG}

  \item{TnBG.ratio}{\link{vector}(optional): a two-element vector containing value and\cr
standard error of ratio of initial Tn signal to BG} 

  \item{rseTn}{\link{numeric}(optional): relative standard error of Tn in percent}

  \item{FR}{\link{vector}(optional): a two-element vector containing value and\cr
standard error of fast ratio of Tn}  

  \item{LnTn.curve}{\link{list}(optional): decay curve data for Ln and Tn, it should contain four elements,\cr 
i.e., \code{names(LnTn.curve)=c("Ln.x","Ln.y","Tn.x","Tn.y")}}

  \item{TxTn}{\link{vector}(optional): ratios of Tx to Tn for various SAR cycles}
}

\details{
Function \link{calED} is used for calculating an equivalent dose and assessing its standard error.
The standard errors of an equivalent dose can be assessd using the Simple Transformation or Monte Carlo method (Duller, 2007).\cr\cr
Interpolation is performed using a combination of golden section search and successive parabolic interpolation (\bold{\emph{R}} function \link{optimize} in package \bold{\emph{stats}}) (freely available Fortran 77 source code at \url{http://www.netlib.org/fmm/fmin.f}). See function \link{fitGrowth} for more details on growth curve fitting.
}

\value{
Return an invisible \link{list} that contains the following elements:
  \item{message}{return 0 if calculation succeeds, 1 if growth curve fitting fails, 2 if natural-dose signal saturates, 
3 if equivalent dose calculation fails, 4 if equivalent dose error assessment fails}
  \item{fitIDX}{Indices of dose points used in growth curve fitting}
  \item{LMpars}{optimized parameters for the growth curve}
  \item{value}{minimized objective for the growth curve}
  \item{avg.error}{average fit error for the growth curve}
  \item{RCS}{reduced chi-square value for the growth curve}
  \item{FOM}{figure of merit value for the growth curve in percent}
  \item{calED.method}{method used for equivalent dose calculation, i.e.,\cr
\code{"Interpolation"} or \code{"Extrapolation"}}
  \item{mcED}{randomly simulated equivalent doses}
  \item{ED}{calculated equivalent dose and its standard error}
  \item{ConfInt}{68 percent and 95 percent confidence intervals for the equivalent dose}
  \item{RecyclingRatio1}{the first recycling ratio and its standard error}
  \item{RecyclingRatio2}{the second recycling ratio and its standard error}
  \item{RecyclingRatio3}{the third recycling ratio and its standard error}
  \item{Recuperation1}{the first recuperation (i.e., ratio of the sensitivity-corrected zero-dose signal to natural-dose signal) 
and its standard error in percent}
  \item{Recuperation2}{the second recuperation (i.e., ratio of the sensitivity-corrected zero-dose signal to maximum regenerative-dose signal) 
and its standard error in percent}
}

\references{
Duller GAT, 2007. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. Ancient TL, 25(1): 15-24.

Duller GAT, 2016. Analyst (v4.31.9), User Mannual.

Galbraith RF, Roberts RG, 2012. Statistical aspects of equivalent dose and error calculation and display in OSL dating: an overview and some recommendations. Quaternary Geochronology, 11: 1-27.
}

\seealso{
  \link{analyseBINdata}; \link{fitGrowth}; \link{calRcyRcp};\cr
  \link{calSARED}; \link{fastED}; \link{calSGCED}
}


\note{
Arguments \code{agID}, \code{Tn}, \code{Tn3BG}, \code{TnBG.ratio}, \code{rseTn}, \code{FR}, \code{LnTn.curve}, 
and \code{TxTn} have nothing to do with equivalent dose calculation. They are used only for plotting purpose.\cr\cr
Argument \code{Tn3BG} indicates if Tn (after background subtraction) is more than 3 sigma above BG, while argument \code{TnBG.ratio} denotes the ratio of Tn (no background subtraction) to BG.\cr\cr
Function \link{calED} will return \code{message=3} (i.e.,\code{"Failed in equivalent dose calculation"}) if the equivalent dose to be calculated below -50 (<Gy>|<s>).\cr\cr
68 percent (one sigma) and 95 percent (two sigma) confidence intervals of equivalent doses will be determined respectively using normal approximation and Monte Carlo simulation,\cr
for \code{errMethod="sp"} and \code{errMethod="mc"}.\cr\cr
Function \code{sgcED} in previous versions was bundled to function \link{calSGCED}.
}

\examples{
  ### Example 1:
  Curvedata <- cbind(c(0, 18, 36, 54, 72, 0, 18),               
                    c(0.026, 1.55, 2.39, 3.46, 4.13, 0.023, 1.61),  
                    c(0.005, 0.11, 0.27, 0.22, 0.20, 0.008, 0.24))                         
  Ltx <- c(3.1,0.31)
  calED(Curvedata, Ltx, model="exp", origin=FALSE)
  
  ### Example 2 (not run):
  # data(BIN)
  # obj_pickBIN <- pickBINdata(BIN, Position=48, 
  #                            LType="OSL", view=FALSE)
  # obj_analyseBIN <- analyseBINdata(obj_pickBIN, nfchn=3, nlchn=20)
  # Curvedata <- obj_analyseBIN$SARdata[-1,3:5]
  # Ltx <- as.numeric(obj_analyseBIN$SARdata[1,4:5])
  # calED(Curvedata, Ltx, model="gok", origin=FALSE)
 
}
\keyword{growth curve}
\keyword{equivalent dose}
