\name{fitGrowth}
\alias{fitGrowth}
\alias{fitGrowth.default}
\title{
Growth curve fitting
}
\description{
Fitting growth curves using the Levenberg-Marquardt algorithm.
}
\usage{
fitGrowth(Curvedata, model = "gok", origin = FALSE, 
          weight = TRUE, trial = FALSE, plot = TRUE, 
          agID = NULL)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): a three-column matrix (i.e., regenerative doses,\cr
sensitivity-corrected regenerative-dose signals, and associated standard errors)}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see details}

  \item{origin}{\link{logical}(optional): logical value indicating if the growth curve should be forced to pass the origin}

  \item{weight}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using a weighted procedure, see details}

  \item{trial}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using other models if the given model fails, see details}

  \item{plot}{\link{logical}(with default): logical value indicating if the results should be plotted}

  \item{agID}{\link{vector}(optional): a three-elemenet vector indicating aliquot (grain) ID, i.e.,\cr
\code{agID[1]=NO}, \code{agID[2]=Position}, \code{agID[3]=Grain}}
}

\details{
In growth curve fitting using function \link{fitGrowth}, five models are available: \cr
\bold{(1)} \emph{"line"}: a linear model, \code{y=a*x+b};\cr
\bold{(2)} \emph{"exp"}: a single saturation exponential model, \code{y=a*[1-exp(-b*x)]+c};\cr
\bold{(3)} \emph{"lexp"}: a single saturation exponential plus linear model, \code{y=a*[1-exp(-b*x)]+c*x+d};\cr
\bold{(4)} \emph{"dexp"}: a double saturation exponential model, \code{y=a*[1-exp(-b*x)]+c*[1-exp(-d*x)]+e};\cr
\bold{(5)} \emph{"gok"}: a general order kinetic model (Guralnik et al., 2015), \code{y=a*[1-(1+b*c*x)^(-1/c)]+d}.\cr\cr
The fitting process is performed using the Levenberg-Marquardt algorithm (minpack: Fortran 90 source code by John Burkardt, freely available at \url{http://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html}). If \code{weight=TRUE}, a weighted procedure will be performed through weighting each data point by its inverse variance. User is advised to set argument \code{plot=TRUE} if possible to visualize the quality of fit.\cr\cr
Argument \code{trial=TRUE} means that if the growth curve can not be fitted successfully using the user-supplied model, then the procedure will try to fit other models instead:
\tabular{lll}{
\bold{Model}              \tab  \bold{Tried model}                       \cr
\code{"dexp"}             \tab  \code{c("dexp", "gok", "exp", "line")}   \cr
\code{"lexp"}             \tab  \code{c("lexp", "gok", "exp", "line")}   \cr
\code{"gok"}              \tab  \code{c("gok", "exp", "line")}           \cr
\code{"exp"}              \tab  \code{c("exp", "line")}                  \cr 
\code{"line"}             \tab  \code{c("line")}                                         
}
For example, if \code{model="dexp"} and \code{trial=TRUE}, then a number of models from sequence\cr
\code{c("dexp", "gok", "exp", "line")} will be applied one after another until the fit succeeds.\cr

The required number of data points for each model is (value inside the parentheses denotes the required number of observations if the model passes the origin):
\tabular{lll}{
\bold{Model}   \tab \bold{Required NPoints}   \cr
\code{"dexp"}  \tab   >=5 (>=4)               \cr
\code{"lexp"}  \tab   >=4 (>=3)               \cr
\code{"gok"}   \tab   >=4 (>=3)               \cr  
\code{"exp"}   \tab   >=3 (>=2)               \cr   
\code{"line"}  \tab   >=2 (>=1)                                                   
}
If user-provided data is not enough for model fitting (i.e., the number of data points is less than the number of parameters of a given model), 
then a model with less parameters will be fitted. For example, if 4 data points are fitted using a \code{"dexp"} (\code{origin=FALSE}) model that actually needs at least 5 data points, 
then a 4-parameter \code{"gok"} model will be fitted instead. 
}

\value{
Return a \link{list} that contains the following elements:
  \item{message}{return 0 if fit succeeds, else 1}
  \item{LMpars}{optimized parameters for the growth curve}
  \item{value}{minimized objective for the growth curve}
  \item{avg.error}{average fit error for the growth curve}
  \item{RCS}{reduced chi-square value for the growth curve}
  \item{FOM}{figure of merit value for the growth curve in percent}
}

\references{
Balian HG, Eddy NW, 1977. Figure-of-merit (FOM), an improved criterion over the normalized chi-squared test for assessing 
goodness-of-fit of gamma-ray spectral peaks. Nuclear Instruments and Methods, 145(2): 389-95.

Guralnik B, Li B, Jain M, Chen R, Paris RB, Murray AS, Li SH, Pagonis V, Valla PG, Herman F, 2015. Radiation-induced growth and isothermal decay of 
infrared-stimulated luminescence from feldspar. Radiation Measurements, 81: 224-231. 

More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory" in Lecture Notes in Mathematics: Numerical Analysis, 
Springer-Verlag: Berlin. 105-116.
}

\note{
Argument \code{agID} has nothing to do with growth curve fitting.\cr\cr
The model to be optimized should not be underdetermined. This means that the number of data points should exceed (or at least be equal to) the number of parameters. For a given model, the procedure will return an error if any standard errors of parameters cannot be estimated by numerical difference-approximation.\cr\cr
}

\seealso{
  \link{analyseBINdata}; \link{SARdata};\cr
  \link{calED}; \link{calSARED}; \link{fastED};\cr
  \link{pickSARdata}; \link{lsNORM}; \link{calSGCED}
}

\examples{
 ### Example 1:
 Curvedata<-cbind(c(0, 18, 36, 54, 72, 0, 18),               
                  c(0.026, 1.55, 2.39, 3.46, 4.13, 0.023, 1.61),  
                  c(0.005, 0.11, 0.27, 0.22, 0.20, 0.008, 0.24)) 
 fitGrowth(Curvedata, model="gok", origin=FALSE, plot=TRUE)

 ### Example 2 (not run):
 # data(SARdata)
 # Curvedata <- SARdata[!is.element(SARdata[,2], "N"),3:5]
 # fitGrowth(Curvedata, model="gok", origin=FALSE)
 
}

\keyword{growth curve}
