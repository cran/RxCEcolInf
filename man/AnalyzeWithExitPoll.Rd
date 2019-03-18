\name{AnalyzeWithExitPoll}
\alias{AnalyzeWithExitPoll}
\title{Workhorse Function for Ecological Inference for Sets of R x C
  Contingency Tables When Incorporating a Survey such as an Exit Poll}
\description{
  This function (using the tuned parameters from \code{TuneWithExitPoll}) fits a
  hierarchical model to data from two sources:  (i) ecological data and in which the underlying
  contigency tables can have any number of rows or columns, and (ii)
  data from a survey sample of some of the contingency tables.  The user
  supplies the data and may specify hyperprior values.  Samples from the
  posterior distribution are returned as an \code{mcmc} object, which can be
  analyzed with functions in the \code{coda} package.
}

\usage{
AnalyzeWithExitPoll(fstring, rho.vec, exitpoll, data = NULL, 
	            num.iters = 1e+06, save.every = 1000, burnin = 10000, 
                    mu.vec.0 = rep(log((0.45/(mu.dim - 1))/0.55), mu.dim), 
                    kappa = 10, nu = (mu.dim + 6), psi = mu.dim, 
                    mu.vec.cu = runif(mu.dim, -3, 0), NNs.start = NULL,  
                    MMs.start = NULL, THETAS.start = NULL, sr.probs = NULL, 
                    sr.reps = NULL, keep.restart.info = FALSE, 
                    keepNNinternals = 0, keepTHETAS = 0, nolocalmode = 50, 
                    numscans = 1, Diri = 100, dof = 4, eschew = FALSE, 
                    print.every = 10000, debug = 1)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{fstring}{String: model formula of contingency tables' column
  totals versus row totals. Must be in specified format (an R
  character string and NOT a true R formula). See Details and
  Examples.}

  \item{rho.vec}{Vector of dimension \eqn{I} = number of contigency
    tables = number of rows in \code{data}: multipliers (usually in
    (0,1)) to the covariance matrix of the proposal distribution for
    the draws of the intermediate level parameters. Typically set to
    the vector output from \code{TuneWithExitPoll}.}
  

  \item{exitpoll}{Matrix of dimensions \eqn{I} = number of contingency
    tables = number of rows in \code{data} by (R * C) = number of cells
    in each contingency table:  The results of a survey sample of
    some of the contingency tables.  Must be in specified format.  See Details.}

  \item{data}{Data frame.}
  
  \item{num.iters}{Positive integer: The number of MCMC iterations for
    the sampler.}
  
  \item{save.every}{Positive integer: The interval at which the draws
    will be saved. \code{num.iters} must be divisible by this value.
    Akin to \code{thin} in some packages. For example, \code{num.iters
    = 1000} and \code{save.every = 10} outputs every 10th draw for a
    total of 100 outputed draws.}
  
  \item{burnin}{Positive integer: The number of burn-in iterations for
  the sampler.}
  
  \item{mu.vec.0}{Vector:  mean of the (normal) hyperprior distribution for
    the \eqn{\mu}{mu} parameter.}
  
  \item{kappa}{Scalar:  The diagonal of the covariance matrix for the
    (normal) hyperprior distribution for the \eqn{\mu}{mu} parameter.}
  
  \item{nu}{Scalar:  The degrees of freedom for the (Inverse-Wishart)
    hyperprior distriution for the \code{SIGMA} parameter.}
  
  \item{psi}{Scalar:  The diagnoal of the matrix parameter of the
    (Inverse-Wishart) hyperprior distribution for the \code{SIGMA} parameter.}
  
  \item{mu.vec.cu}{Vector of dimension \eqn{R*(C-1)}, where
    \eqn{R}(\eqn{C}) is the number of rows(columns) in each contigency
    table:  Optional starting values for \eqn{\mu}{mu} parameter.}
  
  \item{NNs.start}{Matrix of dimension \eqn{I} x (\eqn{R*C}),
    where \eqn{I} is the number of contingency tables = number of rows in
    \code{data}:  Optional starting values for the internal cell
    counts, which must total to the continency table row and column
    totals contained in \code{data}.  Use of the default (randomly
    generated internally) recommended.}
  
  \item{MMs.start}{Matrix of dimension \eqn{I} x (\eqn{R*C}),
    where \eqn{I} is the number of contingency tables = number of rows in
    \code{data}:  Optional starting values for the \emph{missing} internal cell
    counts, which must total to the continency table row and column
    totals contained in \code{data}.  By \emph{missing} internal cell
    counts we mean the counts in the contigency tables not observed in
    the survey sample or \code{exitpoll}.  Use of the default (randomly
    generated internally) recommended.}

  \item{THETAS.start}{Matrix of dimension \eqn{I} x (\eqn{R*C}),
    where \eqn{I} is the number of contingency tables = number of rows in
    \code{data}:  Optional starting values for the contingency table row
    probability vectors.  The elements in each row of \code{THETAS.start}
    must meet \eqn{R} sum-to-one constraints.  Use of the
    default (randomly generated internally) recommended.}
  
  \item{sr.probs}{Matrix of dimension \eqn{I} x \eqn{R}:  Each value
    represents the probability of selecting a particular
    contingency table's row as the row to be calculated deterministically
    in (product multinomial) proposals for Metropolis draws of the
    internal cell counts.  For example, if R = 3 and row 2 of position
    \code{sr.probs} = c(.1, .5, .4), then in the third contingency table
    (correspoding to the third row of \code{data}), the proposal
    algorithm for the interior cell counts will calculate the third
    contingency table's first row deterministically with probability
    .1, the second row with probability .5, and the third row with
    probability .4.  Use of default (generated
    internally) recommended.}
  
  \item{sr.reps}{Matrix of dimension \eqn{I} x \eqn{R}:  Each value
    represents the number of times the (product multinomial proposal)
    Metropolis algorithm will be attempted when, in drawing the internal
    cell counts, the proposal for the
    corresponding contingency table row is to be calculated
    deterministically.  sr.reps has the same structure as sr.probs,
    \emph{i.e.}, position [3,1] of sr.reps corresponds to the third
    contingency table's first row.  Use of default (generated internally) recommended.}
  
  \item{keep.restart.info}{Logical:  Whether last state of the chain
    should be saved to allow restart in the same state.  Restart
    function not currently implemented.}
  
  \item{keepNNinternals}{Positive integer:  The number of draws of the
    internal cell counts in the contingency tables to be outputted.
    Must be divisible into \code{num.iters}.  Use with caution:  results
  in large RAM use even in modest-sized datasets.}
  
  \item{keepTHETAS}{Positive integer:  The number of draws of the
    contingency table
    row probability vectors in the contingency tables to be outputted.
    Must be divisible into \code{num.iters}.   Use with caution:  results
    in large RAM use even in modest-sized datasets.}
  
  \item{nolocalmode}{Positive integer:  How often an alternative drawing
  method for the contigency table internal cell counts will be used.
  Use of default value recommended.}
  
  \item{numscans}{Positive integer:  How often the algorithm to draw the
  contingency table internal cell counts will be implemented before new
  values of the other parameters are drawn.  Use of default value recommended.}
  
  \item{Diri}{Positive integer:  How often a product Dirichlet proposal
    distribution will be used to draw the contingency table row
    probability vectors (the THETAS).}
  
  \item{dof}{Positive integer:  The degrees of freedom of the
    multivariate \eqn{t} proposal distribution used in drawing the contingency
    table row probability vectors (the THETAS).}
  
  \item{eschew}{Logical:  If true, calculation of certain functions of
    the interntal cell counts omits the two right-most columns instead
    of only the single right-most column.  Not yet implemented.}
  
  \item{print.every}{Positive integer:  If \code{debug == 1}, the number
    of every \code{print.every}\eqn{th} iteration will be written to the
    screen.  Must be divisible into \code{num.iters}.}
  
  \item{debug}{Integer:  Akin to \code{verbose} in some packages.  If set
    to 1, certain status information (including rough notification
    regarding the number of iterations completed) will be
  written to the screen.}
}
\details{
  \code{AnalyzeWithExitPoll} is the workhorse function in fitting the R
  x C ecological inference model described in Greiner & Quinn (2009)
  with the addition of information from a survey sample from some of the
  contingency tables.  Details and terminology of the basic
  (\emph{i.e.}, without a survey sample) data structure and ecological
  inference model are discussed in the
  documentation accompanying the \code{Analyze} function.  

  In the present implementation, the \code{AnalyzeWithExitPoll}
  presumes that the survey consisted of a simple random sample from the
  in-sample contingency tables.  Future implementations will allow
  incorporation of more complicated survey sampling schemes.

  The arguments to \code{AnalyzeWithExitPoll} are essentially identical
  to those of \code{Analyze} with the major exception of
  \code{exitpoll}.  \code{exitpoll} feeds the results of the survey
  sample to the function, and a particular format is required.
  Specifically, \code{exitpoll} must have the same number of rows as
  \code{data}, meaning one row for each contigency table in the
  dataset.  It must have R * C columns, meaning one column for each cell
  in one of the ecological data's contingency tables.  The first row of
  \code{exitpoll} must correspond to the first row of \code{data},
  meaning that the two rows must contain information from the same
  contingency table.  The second row of \code{exitpoll} must contain
  information from the contingency table represented in the second row
  of \code{data}.  And so on.  Finally, \code{exitpoll} must have counts
  from the sample of the contingency table in vectorized row major
  format.

  To illustrate with a voting example:  Suppose the contingency tables have two rows, labeled
  bla and whi, and three columns, denoted Dem, Rep, and Abs.
  In other words, the \code{fstring} argument would be \code{"Dem, Rep, Abs ~
  bla, whi"}.  Suppose there are 100 contingency tables.  The \code{data}
  will be of dimension \eqn{100 \times 5}, with each row consisting of
  the row and column totals from that particular contigency table.
  \code{exitpoll} will be of dimension \eqn{100 \times 6}.  Row 11 of
  the \code{exitpoll} will consist of the following:  in position 1, the
  number of blacks voting Democrat observed in the sample of contingency
  table 11; in position 2, the number of blacks voting Republican
  observed in the sample of contigency table 11; in position 3, the
  number of blacks Abstaining from voting observed in the sample of
  contingency table 11; in position 4, the number of whites voting
  Democrat observed in the sample of contingency table 11; etc.

  For tables in which there was no sample taken (\emph{i.e.},
  out-of-sample tables), the corresponding row of \code{exitpoll} should have a vector of 0s.

  Model fitting proceeds similarly as in \code{Analyze}, and output is
  simimilarly similar.  See documentation accompanyng \code{Analyze} for
  further information.
}

\value{
  An object of class \code{mcmc} suitable for use in functions in the coda
  package.  Additional items, listed below, may be retrieved from this object, as
  detailed in the examples section.

  \item{dim}{Vector (integers) of length 2: number of saved
  simulations and number of automatically outputted parameters.}

  \item{dimnames}{List: the first element \code{NULL} (currently not
  used), and the second element is a vector of the names of
  the automatically outputted parameters.}

  \item{acc.t}{Vector of length \eqn{I} = number of contigency tables:  The
    fraction of multivariate \eqn{t} proposals accepted in the Metropolis
    algorithm used to draw the \code{THETA}s (meaning the intermediate
    parameters in the hierarchy).}

  \item{acc.Diri}{Vector of length \eqn{I} = number of contigency tables:  The
    fraction of Dirichlet-based proposals accepted in the Metropolis
    algorithm used to draw the \code{THETA}s (meaning the intermediate
    parameters in the hierarchy).}

  \item{vld.multinom}{Matrix:  To draw from the conditional posterior of
    the internal cell counts of a contigency table, the \code{Analyze} function
    draws R-1 vectors of lenth C from multinomial distributions.  In
    then calculates the counts in the additional row (denote this row as
    r') deterministically.  This procedure can result in negative values
    in row r', in which case the overall proposal for the interior cell
    counts is outside the parameter space (and thus invalid).
    vld.multinom keeps track of the percentage of proposals drawn in
    this manner that are valid (\emph{i.e.}, not invalid).  Each row of
    vld.multinom corresponds to a
    contingency table.  Each column in vld.multinom
    corresponds to a \emph{row} in the a contingency table.  Each entry
    specifies the percentage of multinomial proposals that are valid
    when the specified contingency table row serves as the r' row.  For
    instance, in position 5,2 of vld.multinom is the fraction of valid
    proposals for the 5th contingency table when the second contigency
    table row is the r'th row.  A value of ``NaN'' means that \code{Analyze}
    chose to use a different (slower) method of drawing the internal
    cell counts because it suspected that the multinomial method would
    behave badly.}

  \item{acc.multinom}{Matrix:  Same as vld.multinom, except the entries
    represent the fraction of proposals accepted (instead of the
    fraction that are in the permissible parameter space).}

  \item{numrows.pt}{Integer:  Number of rows in each contingency table.}

  \item{numcols.pt}{Integer:  Number of columns in each contingency
    table.}

  \item{THETA}{mcmc:  Draws of the \code{THETA}s.  See Details and Examples.}

  \item{NN.internals}{mcmc:  Draws of the internal cell counts.  See
    Details and Examples.}

}
\references{
  D. James Greiner \& Kevin M. Quinn. 2009.  ``R x C Ecological
  Inference:  Bounds, Correlations, Flexibility, and Transparency of
  Assumptions.'' \emph{J.R. Statist. Soc. A} 172:67-81.

  Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
}
\author{D. James Greiner, Paul D. Baines, \& Kevin M. Quinn}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\section{Warnings}{ \strong{Computer time:} At present, using this
  function (and the others in this package) requires substantial
  computer time. The lack of information in ecological data results in
  slow mixing chains, and the number of parameters that must be drawn
  in each Gibbs sampler iteration is large. Chain length should be
  adjusted to achieve adequate convergence. In general, the more
  segregated the housing patterns in the jurisdiction (meaning the
  greater the percentage of contingency tables in which one row's
  counts make up a large portion of that table's total), the smaller
  the number of iterations needed. We are exploring more efficient
  sampling algorithms that we anticipate will result in better mixing
  and faster drawing. At present, however, users should anticipate
  that analysis of a dataset will take several hours.
  
  \strong{Large datasets:} At present, use of this fuction (and thus
    this package) is not recommended for large (\emph{i.e.}, more than
    1000 contingency tables) datasets. See immediately above.
  
  \strong{RAM requirements:} Do not select large values of
    \code{keepNNinternals} or \code{keepTHETAS} without adequate RAM.
  
  \strong{Gelman-Rubin diagnostic in the coda package:} Using the
    Gelman-Rubin convergence diagnostic as presently implemented in
    the CODA package (called by \code{gelman.diag()}) on multiple
    chains produced by Analyze will cause an error. The reason is that
    some of the NN.internals and functions of them
    (\eqn{LAMBDAs}{\Lambda's}, \code{TURNOUTs},
    \eqn{GAMMAs}{\Gamma's}, and \eqn{BETAs}{\beta's}) are linearly
    dependant, and the current coda implmentation of gelman.diag()
    relies on a matrix inversion.}


%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
\dontrun{
SimData <- gendata.ep()    #  simulated data
FormulaString <- "Dem, Rep, Abs ~ bla, whi, his"
EPInvTune <-  TuneWithExitPoll(fstring = FormulaString,
                               data = SimData$GQdata,
                               exitpoll=SimData$EPInv$returnmat.ep,
                               num.iters = 10000,
                               num.runs = 15)
EPInvChain1 <- AnalyzeWithExitPoll(fstring = FormulaString,
                                   data = SimData$GQdata,
                                   exitpoll=SimData$EPInv$returnmat.ep,
                                   num.iters = 2000000,
                                   burnin = 200000,
                                   save.every = 2000,
                                   rho.vec = EPInvTune$rhos,
                                   print.every = 20000,
                                   debug = 1,
                                   keepTHETAS = 0,
                                   keepNNinternals = 0)
EPInvChain2 <- AnalyzeWithExitPoll(fstring = FormulaString,
                                   data = SimData$GQdata,
                                   exitpoll=SimData$EPInv$returnmat.ep,
                                   num.iters = 2000000,
                                   burnin = 200000,
                                   save.every = 2000,
                                   rho.vec = EPInvTune$rhos,
                                   print.every = 20000,
                                   debug = 1,
                                   keepTHETAS = 0,
                                   keepNNinternals = 0)
EPInvChain3 <- AnalyzeWithExitPoll(fstring = FormulaString,
                                   data = SimData$GQdata,
                                   exitpoll=SimData$EPInv$returnmat.ep,
                                   num.iters = 2000000,
                                   burnin = 200000,
                                   save.every = 2000,
                                   rho.vec = EPInvTune$rhos,
                                   print.every = 20000,
                                   debug = 1,
                                   keepTHETAS = 0,
                                   keepNNinternals = 0)
EPInv <- mcmc.list(EPInvChain1, EPInvChain2, EPInvChain3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
