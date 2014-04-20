\name{calED}
\alias{calED}
\alias{calED.default}
\title{
Fitting a dose response curve and calculating some equivalent dose values
}
\description{
Fitting a dose response curve and calculating some equivalent dose values, the error of an equivalent dose is assessed with methods described in Duller GAT (2007a).
}
\usage{
calED(Curvedata, Ltx, 
      model = c("line","exp","line+exp"), origin = FALSE,
      nstart = 100, upb = 1, ErrorMethod = c("mc", "sp"),
      nsim = 1000, plot = TRUE, outfile = NULL)
}
\arguments{
  \item{Curvedata}{\link{data.frame}(\bold{required}): three columns, i.e. regenerative doses (redose1, redose2, etc), standardised OSL signals (Lx1/Tx1, Lx2/Tx2, etc), and standard errors of 
standardised OSL signals}
  \item{Ltx}{\link{vector} or \link{matrix}(\bold{required}): a vector contains two values (corresponding to the standardised OSL signal and its standard error) or a matrix contains two 
columns (each row corresponds to a standardised OSL signal and an associated standard error in an application of the SGC method)}
  \item{model}{\link{character}(with default): a model used for fitting the growth curve ("line", "exp" or "line+exp"), default model="line". See details}
  \item{origin}{\link{logical}(optional): force the fitting to pass the origin (x=0,y=0) or not}
  \item{nstart}{\link{numeric}(with default): the maximum number of attempts of initializations in the non-linear growth curve fitting}
  \item{upb}{\link{numeric}(with default): the upper limit on b value, initial b value is generated uniformly from the space (0,upb), see details}
  \item{ErrorMethod}{\link{character}(with default): method ("sp" or "mc") used for assessing the standard error of an equivalent dose. Default ErrorMethod="mc", see details}
  \item{nsim}{\link{numeric}(with default): the allowed maximum number of simulations for applying a Monte Carlo method to assess the standard error of an equivalent dose}
  \item{plot}{\link{logical}(with default): draw the fitted dose response curve or not}
  \item{outfile}{\link{character}(with default): if specified, simulated ED values will be written to a file named "outfile" in ".csv" format and will be saved to the current work directory}
}
\details{
This function can be used for both dose response curve fitting and equivalent dose calculating. Calculating a single ED value or a number of ED values (SGC method) are possible, three fitting methods are available: \cr
1) "line": a linear model, y=a*x+b \bold{or} y=a*x if origin=TRUE;\cr
2) "exp": an exponential model, y=a*(1-exp(-b*x))+c \bold{or} y=a*(1-exp(-b*x)) if origin=TRUE;\cr
3) "line+exp": an exponential plus linear model, y=a*(1-exp(-b*x))+c*x+d \bold{or} y=a*(1-exp(-b*x))+c*x if origin=TRUE.\cr\cr
Standard error of an equivalent dose can be assessd using two method (Duller GAT, 2007a): \cr
a) "sp": a simple transformation of s(Lx/Tx) to s(ED); \cr
b) "mc": a Monte Carlo method \cr\cr
Curve fitting is performed using the Levenberg-Marquardt algorithm (minpack: original FORTRAN 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. FORTRAN 90 version by John Burkardt). Interpolation is performed using a combination of golden section search and successive parabolic interpolation (the same as function \link{optimize} in package "stats").\cr\cr
As for the fitting of a linear model, initialization is not a problem, but if the model to be fitted is either an exponential or an exponential plus linear model, then a series of b values are initialized randomly with a uniform distribution in the space (0,upb), then other parameters (a or a, c or a, c, d) can be estimated with a Linear Algebra method. User is advised to set the default setting "plot" to be TRUE if possible in order to see if the fitted growth curve is correct.
}
\value{
Return an invisible list that contains following elements:
  \item{LMpars}{the estimated characteristic parameters and standard errors of the growth curve}
  \item{value}{the sum of squared residuals}
  \item{fit.value}{observations and fitted values}
  \item{ED}{the estimated equivalent dose (and its standard error), which may be a two-element vector or a two-column matrix, depending on argument "Ltx"}
}
\references{
Roberts, H.M. and Duller, G.A.T., 2004. Standardised growth curves for optical dating of sediment using multiple-grain aliquots. Radiation Measurements 38(2), pp. 241-252.

Duller, G.A.T., 2007a. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. 
Ancient TL 25(1), pp. 15-24.

Duller, G., 2007b. Analyst. pp. 27-28.
}
\note{
Interpolation is performed within the maximum regenerative dose, so the value of Ltx that used to calculate an equivalent dose cannot be larger than the maximum standardised signals that used to construct the growth curve if the growth curve is fitted with a non-linear model. Also note that if any estimated characteristic parameter of the growth curve is below or equal to zero, the estimates will be suppressed and will be regarded as a fail.
}
\seealso{
  \link{optimize}
}
\examples{
  # load package "numOSL"
  library(numOSL)
  # Example 1: calculating equivalent doses using data from Duller (2007a)
  data(Curvedata)
  Res<-matrix(nrow=4,ncol=3)
  for (i in 1:4)  {
    model<-ifelse(i==1 || i==2,"exp","line")    
    Simple<-calED(Curvedata[[i]][[1]],Curvedata[[i]][[2]],
                  model=model,ErrorMethod="sp",plot=FALSE)
    MonteCarlo<-calED(Curvedata[[i]][[1]],Curvedata[[i]][[2]],
                      model=model,ErrorMethod="mc",plot=FALSE)
    Res[i,1:2]<-Simple$ED
    Res[i,3]<-MonteCarlo$ED[2]
  } # end for
  colnames(Res)<-c("ED","Simp.Error","MC.Error")
  rownames(Res)<-c("Exponential","Exponential","linear","linear")
  Res
  # Example 2: performing a SGC method
  ltx1<-c(0.032,1.61,2.55,3.21,3.87,0.031,1.55) # Lx/Tx for the first aliquot
  ltx2<-c(0.025,1.44,2.47,3.35,4.17,0.033,1.47) # Lx/Tx for the second aliquot
  ltx3<-c(0.027,1.51,2.68,3.52,4.41,0.021,1.39) # Lx/Tx for the third aliquot
  ltx4<-c(0.018,1.71,2.28,3.81,4.03,0.017,1.62) # Lx/Tx for the four aliquot
  ltx5<-c(0.026,1.49,1.99,3.43,4.17,0.015,2.01) # Lx/Tx for the five aliquot
  ltx<-cbind(ltx1,ltx2,ltx3,ltx4,ltx5)
  ltx<-cbind(apply(ltx,MARGIN=1,mean),
             apply(ltx,MARGIN=1,sd))  # means and standard deviations
  redose<-c(0,18,36,54,72,0,18)       # the same ReDose for the five aliquots
  Curvedata<-data.frame(redose,ltx)
  Ltx<-cbind(c(0.5,1.0,1.8,2.3,2.8,3.1,3.6,4.0),rep(0.1,8))         
  calED(Curvedata,Ltx,model="exp",origin=FALSE) # fitting y=a*(1-exp(-b*x))+c
}
\keyword{equivalent dose}
