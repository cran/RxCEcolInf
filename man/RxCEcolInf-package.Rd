\name{RxCEcolInf-package}
\alias{RxCEcolInf-package}
\alias{RxCEcolInf}
\docType{package}
\title{
 RxCEcolInf
}
\description{
Fits the R x C ecological inference model described in Greiner and Quinn
(2009).  Allows the inclusion of survey information.
}
\details{
\tabular{ll}{
Package: \tab RxCEcolInf\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2009-04-14\cr
License: \tab \cr
}
The user should place the data from the contingency tables into a
dataframe, each line of which represents one table.  The function Tune
should be called first; this will tune the markov chain monte carlo
algorithm used to fit the model.  The user feeds the results from Tune
into Analyze, which produces the results.  Both Tune and Analyze are
called using a string that resembles the R formula interface.  If a
simple random sample is available from certain contingency tables, the
user may incorporate this survey using TuneWithExitPoll and AnalyzeWithExitPoll.
}
\author{
D. James Greiner, Paul D. Baines, \& Kevin M. Quinn

Maintainer: Kevin M. Quinn <kquinn@fas.harvard.edu>
}
\references{
   D. James Greiner \& Kevin M. Quinn. 2009.  ``R x C Ecological
  Inference:  Bounds, Correlations, Flexibility, and Transparency of
  Assumptions.'' \emph{J.R. Statist. Soc. A} 172:67-81.
}
\keyword{ models }
\examples{
\dontrun{
library(RxCEcolInf)
data(stlouis)
Tune.stlouis <- Tune("Bosley, Roberts, Ribaudo, Villa, NoVote ~ bvap, ovap",
                     data = stlouis,
                     num.iters = 10000,
                     num.runs = 15)
Chain1.stlouis <- Analyze("Bosley, Roberts , Ribaudo, 
                          Villa, NoVote ~ bvap, ovap",
                          rho.vec = Tune.stlouis$rhos,
                          data = stlouis,
                          num.iters = 1500000,
                          burnin = 150000,
                          save.every = 1500,
                          print_every = 15000,
                          debug = 1,
                          keepNNinternals = 100,
                          keepTHETAS = 100)
Chain2.stlouis <- Analyze("Bosley, Roberts , Ribaudo, Villa, 
                          NoVote ~ bvap, ovap",
                          rho.vec = Tune.stlouis$rhos,
                          data = stlouis,
                          num.iters = 1500000,
                          burnin = 150000,
                          save.every = 1500,
                          print_every = 15000,
                          debug = 1,
                          keepNNinternals = 100,
                          keepTHETAS = 100)
Chain3.stlouis <- Analyze("Bosley, Roberts , Ribaudo, Villa, 
	                  NoVote ~ bvap, ovap",
                          rho.vec = Tune.stlouis$rhos,
                          data = stlouis,
                          num.iters = 1500000,
                          burnin = 150000,
                          save.every = 1500,
                          print_every = 15000,
                          debug = 1,
                          keepNNinternals = 100,
                          keepTHETAS = 100)
stlouis.MCMClist <- mcmc.list(Chain1.stlouis, Chain2.stlouis,
Chain3.stlouis)
}
}
