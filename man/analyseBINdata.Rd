\name{analyseBINdata}
\alias{analyseBINdata}
\alias{analyseBINdata.default}
\alias{analyseBINdata0}
\alias{analyseBINdata0.default}

\title{BIN data analysis}

\description{Analysing signal data records extracted from a BIN file.}

\usage{
analyseBINdata(obj_pickBIN, nfchn, nlchn, bg = "late", 
               me = 2.0, distp = "p", kph = NULL, 
               kdc = NULL, dcr = NULL, FR.fchn = NULL, 
               FR.mchn = NULL, FR.lchn = NULL, 
               signal.type = "LxTx", outfile = NULL)

analyseBINdata0(obj_pickBIN, fchn, lchn, bg="late", me=2.0, 
                distp="p", kph=NULL, kdc=NULL, dcr=NULL, 
                FR.fchn=NULL, FR.mchn=NULL, FR.lchn=NULL, 
                signal.type="LxTx", outfile=NULL)
}

\arguments{
  \item{obj_pickBIN}{\link{list}(\bold{required}): an object of S3 class \code{"pickBIN"} produced by\cr 
function \link{pickBINdata}}

  \item{nfchn}{\link{integer}(\bold{required}): number of the first few channels from the initial\cr
part of a decay curve. Number of counts summed over channels\cr
\code{(Delay+1L):(Delay+nfchn)} is calculated as the fast-component\cr
plus backgroud signal}

  \item{nlchn}{\link{integer}(\bold{required}): number of the last few channels from the latter part\cr
of a decay curve. If \code{bg="late"}, number of counts averaged over channels\cr
\code{(Delay+On-nlchn+1L):(Delay+On)} will be calculated as the backgroud\cr
signal, if \code{bg="early"}, number of counts averaged over channels\cr
\code{(Delay+nfchn+1L):(Delay+nfchn+nlchn)} will be calculated as the\cr
backgroud signal. \code{Delay} and \code{On} are obtained internally from the BIN file. 
}

\item{fchn}{\link{integer}(\bold{required}): channels used for calculating the fast-component signals}

\item{lchn}{\link{integer}(\bold{required}): channels used for calculating the background counts}

  \item{bg}{\link{character}(with default): background subtraction method, i.e.,\cr
\code{bg="early"} or \code{bg="late"}}

  \item{me}{\link{numeric}(with default): measurement error of Lx (or Tx) in percent}

  \item{distp}{\link{character}(with default): distribution of photon counts, \code{distp="p"} denotes\cr
Poisson distribution, \code{distp="op"} denotes Over Poisson distribution}

  \item{kph}{\link{numeric}(optional): correction factor for photon counts}

  \item{kdc}{\link{numeric}(optional): correction factor for dark counts}

  \item{dcr}{\link{numeric}(optional): dark count rate}

  \item{FR.fchn}{\link{vector}(optional): fast-component signal channels, note that those channels are extracted
internally from the "ON" channels,\cr
i.e., \code{FR.fchn=((Delay+1):(Delay+On))[FR.fchn]}.\cr
Example: \code{FR.fchn=1:5}}

  \item{FR.mchn}{\link{vector}(optional): medium-component signal channels, note that those channels are extracted 
internally from the "ON" channels,\cr
i.e., \code{FR.mchn=((Delay+1):(Delay+On))[FR.mchn]}.\cr
Example: \code{FR.mchn=31:60}}

  \item{FR.lchn}{\link{vector}(optional): background signal channels, note that those channels are extracted 
internally from the "ON" channels,\cr
i.e., \code{FR.lchn=((Delay+1):(Delay+On))[FR.lchn]}.\cr
Example: \code{FR.lchn=201:250}}

  \item{signal.type}{\link{character}(with default): type of signal, \code{"LxTx"}, \code{"Lx"}, or \code{"Tx"}}

  \item{outfile}{\link{character}(optional): if specified, analysis results (i.e., \code{NO}, \code{Position}, \code{Grain}, 
\code{SAR.Cycle}, \code{Dose}, \code{Init}, \code{BG}, \code{Lx}, \code{seLx}, \code{TInit}, \code{TBG}, \code{Tx}, \code{seTx}, 
\code{LxTx}, \code{seLxTx}) will be written to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\details{
Function \link{analyseBINdata} is used for signal (i.e., Lx, Tx, and Lx/Tx) calculation. 
It provides two protocols for background subtraction (i.e., the early and late background subtraction methods).\cr\cr
Standard error of signals are assessed using two methods: if photon counts are assummed to follow Poisson distributions, 
\bold{Eqn.(3)} of Galbraith (2002) will be applied; if photon counts are over-dispersed, \bold{Eqn.(10)} of Bluszcz et al. (2015) will be applied.\cr\cr
If arguments \code{FR.fchn}, \code{FR.mchn}, and \code{FR.lchn} are provided, fast ratio will be calculated according to Madsen et al. (2009).
}

\value{
Return an invisible \link{list} of S3 class object "analyseBIN" containing the following elements:
  \item{SARdata}{a \link{data.frame} containing calculated SAR data sets}

  \item{criteria}{values used as rejection criteria (0-1 values indicating if Tn is more than 3 sigma above BG or not, 
ratio of initial Tn signal to BG and associated standard error, relative standard error of Tn in percent, fast ratio of Tn and associated standard error), 
\code{NA} is produced if the value can not be calculated. Note that in this function rejection criteria are calculated but not applied}

  \item{Tn}{values of Tn and associated standard errors}

  \item{LnTn.curve}{decay curves for Ln and Tn for different aliquots (grains)}

  \item{TxTn}{ratios of Tx to Tn for various SAR cycles}

  \item{agID}{aliquot or grain ID (i.e., \code{NO}, \code{Position}, and \code{Grain})}

\code{SARdata} is a \link{data.frame} containing the following elements if \code{signal.type="LxTx"}:
\tabular{lll}{
\bold{Element}    \tab  \bold{Description}                             \cr
\code{NO}         \tab  aliquot (grain) number                         \cr
\code{SAR.Cycle}  \tab  SAR cycle (N, R1, R2, R3, ...)                 \cr 
\code{Dose}       \tab  regenerative dose                              \cr
\code{LxTx}       \tab  sensitivity-corrected regenerative-dose signal \cr
\code{seLxTx}     \tab  standard error of LxTx                                  
}

\code{SARdata} contains the following elements if \code{signal.type="Lx"}:
\tabular{lll}{
\bold{Element}    \tab  \bold{Description}                \cr
\code{NO}         \tab  aliquot (grain) number            \cr
\code{SAR.Cycle}  \tab  SAR cycle (N, R1, R2, R3, ...)    \cr 
\code{Dose}       \tab  regenerative dose                 \cr
\code{Lx}         \tab  regenerative-dose signal          \cr
\code{seLx}       \tab  standard error of Lx                                 
}

\code{SARdata} contains the following elements if \code{signal.type="Tx"}:
\tabular{lll}{
\bold{Element}    \tab  \bold{Description}                   \cr
\code{NO}         \tab  aliquot (grain) number               \cr
\code{SAR.Cycle}  \tab  SAR cycle (N, R1, R2, R3, ...)       \cr 
\code{Dose}       \tab  regenerative dose                    \cr
\code{Tx}         \tab  test-dose signal                     \cr
\code{seTx}       \tab  standard error of Tx                                 
}
}

\references{
Ballarini M, Wallinga J, Wintle AG, Bos AJJ, 2007. A modified SAR protocol for optical dating of individual grains from young quartz samples. 
Radiation Measurements, 42(3): 360-369.

Bluszcz A, Adamiec G, Heer AJ, 2015. Estimation of equivalent dose and its uncertainty in the OSL SAR protocol when count numbers do not follow a Poisson distribution. Radiation Measurements, 81: 46-54. 

Cunningham AC, Wallinga J, 2010. Selection of integration time intervals for quartz OSL decay curves. Quaternary Geochronology, 5(6): 657-666

Duller GAT, 2016. Analyst (v4.31.9), User Mannual.

Durcan JA, Duller GAT, 2011. The fast ratio: A rapid measure for testing the dominance of the fast component in the initial OSL signal from quartz. Radiation Measurements, 46(10): 1065-1072.

Galbraith R, 2002. A note on the variance of a backround-corrected OSL count. Ancient TL, 20(2): 49-51.

Madsen AT, Duller GAT, Donnelly JP, Roberts HM, Wintle AG, 2009. A chronology of hurricane landfalls at Little Sippewissett Marsh, Massachusetts, USA, using optical dating. Geomorphology, 109(1-2): 36-45.
}

\note{
Though function \link{analyseBINdata} is originally designed to analyze CW-OSL data sets, IRSL data sets obtained from the SAR protocol can also be analyzed.
}

\seealso{
  \link{loadBINdata}; \link{pickBINdata}; \link{pickSARdata}; \link{calED};\cr
  \link{calSARED}; \link{calSGCED}; \link{fitGrowth}; \link{lsNORM}; \link{BIN}
}

\examples{
   ### Example 1 (not run):
   # obj_loadBIN <- loadBINdata("foo.bin", view=TRUE)
   # obj_pickBIN <- pickBINdata(obj_loadBIN, Position=2, LType="OSL")
   # analyseBINdata(obj_pickBIN, nfchn=3, nlchn=20)

   ### Example 2:
   data(BIN)
   obj_pickBIN <- pickBINdata(BIN, Position=c(2,4,6,8,10), 
                              LType="OSL", view=FALSE)
   obj_analyseBIN <- analyseBINdata(obj_pickBIN, nfchn=4, nlchn=20) 
   #obj_analyseBIN <- analyseBINdata0(obj_pickBIN, fchn=1:4, nlchn=231:250)
   obj_analyseBIN$SARdata
}

\keyword{growth curve}
\keyword{BIN file}
