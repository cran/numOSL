\name{calED}
\alias{calED}
\alias{calED.default}
\title{
Fitting a dose response curve and calculating an equivalent dose
}
\description{
Fitting a dose response curve and calculating an equivalent dose, the error of the equivalent dose is assessed after Duller GAT (2007a).
}
\usage{
calED(Curvedata, Ltx, 
      model = c("line","exp","line+exp"), origin = FALSE,
      nstart = 100, upb = 1, ErrorMethod = c("mc", "sp"),
      nsim = 1000, plot = TRUE, samplename = NULL)
}
\arguments{
  \item{Curvedata}{\link{data.frame}(\bold{required}): three columns, i.e. regenerative doses (redose1, redose2, etc), standardized signals (Lx1/Tx1, Lx2/Tx2, etc), errors of standardized signals}
  \item{Ltx}{\link{vector}(\bold{required}): two values (the standardized natural signal and its standard error), then an equivalent dose is estimated with interpolation}
  \item{model}{\link{character}(with default): a fitting model ("line", "exp" or "line+exp"), default model="line". See details}
  \item{origin}{\link{logical}(optional): force the fitting to pass the origin (x=0,y=0) or not}
  \item{nstart}{\link{numeric}(with default): the maximum number of attempts that used for initializing parameters of curve fitting}
  \item{upb}{\link{numeric}(with default): the upper limit on b value, initial b value is generated uniformly from the space (0,upb), see details}
  \item{ErrorMethod}{\link{character}(with default): method ("sp" or "mc") used for assessing the standard error of equivalent dose. Default ErrorMethod="mc", see details}
  \item{nsim}{\link{numeric}(with default): the allowed maximum number of simulations for application of a Monte Carlo method to assess the standard error of an equivalent dose}
  \item{plot}{\link{logical}(with default): draw the fitted dose response curve or not}
  \item{samplename}{\link{character}(with default): the name of the sample}
}
\details{
This function can be used for both dose response curve fitting and equivalent dose calculating, three fitting methods are available: \cr
1) "line": a linear model, y=a*x+b \bold{or} y=a*x if origin=TRUE;\cr
2) "exp": an exponential model, y=a*(1-exp(-b*x))+c \bold{or} y=a*(1-exp(-b*x)) if origin=TRUE;\cr
3) "line+exp": an exponential plus linear model, y=a*(1-exp(-b*x))+c*x+d \bold{or} y=a*(1-exp(-b*x))+c*x if origin=TRUE.\cr\cr
Standard error of equivalent dose can be accessed using two method (Duller GAT, 2007a): \cr
a) "sp": a simple transformation of s(Lx/Tx) to s(ED); \cr
b) "mc": a Monte Carlo method \cr\cr
Curve fitting is performed using the Levenberg-Marquardt algorithm (minpack: original FORTRAN 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. FORTRAN 90 version by John Burkardt). Interpolation is performed using a combination of golden section search and successive parabolic interpolation (the same as function \link{optimize} in package "stats").\cr\cr
As for the fitting of a linear model, initial parameters is not a problem, but if the model that to be fitted is either an exponential or an exponential plus linear model, then a series of b values are initialized randomly from a uniform distribution in the space (0,upb), then other parameters (a or a, c or a, c, d) can be estimated with a Linear Algebra method. User is advised to keep the default setting "plot" to be TRUE if possible in order to see whether the fitted dose response curve is correct.
}
\value{
Return an invisible list that contains following elements:
  \item{mcED}{an equivalent doses obtained by Monte Carlo simulation, available only if ErrorMethod="mc", or else it will be NULL}
  \item{LMpars}{the characteristic parameters and standard errors for the dose response curve estimated with the Levenberg-Marquardt method}
  \item{residual}{the sum of square of residuals}
  \item{fit.value}{observations and fitted values}
  \item{ED}{estimated equivalent dose (and standard error)}
}
\references{
Duller, G.A.T., 2007a. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. 
Ancient TL 25(1), pp. 15-24.

Duller, G., 2007b. Analyst. pp. 27-28.
}
\note{
Interpolation is performed within the maximum regenerative dose, so the value of Lt/Tx that used to calculate an equivalent dose must not larger than the maximum standardized signals that used to construct the dose response curve. 
}
\seealso{
  \link{sgcED}; \link{optimize}
}
\examples{
  # load package "numOSL"
  library(numOSL)
  # Calculating equivalent doses using data from Duller (2007a)
  data(Curvedata)
  Res<-matrix(nrow=4,ncol=3)
  par(mfrow=c(2,2))
  for (i in 1:4)  {
    model<-ifelse(any(i==c(1,2)),"exp","line")    
    Simple<-calED(Curvedata[[i]][[1]],Curvedata[[i]][[2]],
                  model=model,ErrorMethod="sp",plot=FALSE)
    MonteCarlo<-calED(Curvedata[[i]][[1]],Curvedata[[i]][[2]],
                      model=model,ErrorMethod="mc",plot=TRUE)
    Res[i,1:2]<-Simple$ED
    Res[i,3]<-MonteCarlo$ED[2]
  } # end for
  colnames(Res)<-c("ED","Simp.Error","MC.Error")
  rownames(Res)<-c("Exponential","Exponential","linear","linear")
  print(Res)
  par(mfrow=c(1,1))
}
\keyword{equivalent dose}
