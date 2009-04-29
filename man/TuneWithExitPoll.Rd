\name{TuneWithExitPoll}
\alias{TuneWithExitPoll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tuning Function for Ecological Inference for Sets of R x C
  Contingency Tables When Incorporating a Survey such as an Exit Poll}
\description{
  This function tunes the markov chain monte carlo algorithm used to fit a
  hierarchical model to data from two sources:  (i) ecological data in which the underlying
  contigency tables can have any number of rows or columns, and (ii data
  from a survey samlpe of some of the contingency tables.  The user
  supplies the data and may specify hyperprior values.  The function's
  primary output is a vector of multipliers, called \code{rhos}, used to
  adjust the covariance matrix of the multivariate \eqn{t_4} distribution used
  to propose new values of intermediate-level parameters (denoted
  THETAS).
}

\usage{
TuneWithExitPoll(fstring, exitpoll, data = NULL, num.runs = 12, 
                 num.iters = 10000, rho.vec = rep(0.05, ntables), 
                 kappa = 10, nu = (mu.dim + 6), psi = mu.dim, 
                 mu.vec.0 = rep(log((0.45/(mu.dim - 1))/0.55), mu.dim), 
                 mu.vec.cu = runif(mu.dim, -3, 0), nolocalmode = 50,  
                 sr.probs = NULL, sr.reps = NULL, numscans = 1, 
                 Diri = 100, dof = 4, debug = 1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fstring}{String: model formula of contingency tables' column
  totals versus row totals. Must be in specified format (an R
  character string and NOT a true R formula).}
  
  \item{exitpoll}{Matrix of dimensions \eqn{I} = number of contingency
    tables = number of rows in \code{data} by (R * C) = number of cells
    in each contingency table:  The results of a survey sample of
    some of the contingency tables.  Must be in specified format.  See Details.}
  
  \item{data}{Data frame.}
  
  \item{num.runs}{Positive integer:  The number of runs or times (each of
    \code{num.iters} iterations) the tuning algorthm will be
    implemented.}
  
  \item{num.iters}{Positive integer:  The number of iterations in each
    run of the tuning algorithm.}
  
  \item{rho.vec}{Vector of dimension \eqn{I} = number of contigency
    tables = number of rows in \code{data}: initial values of multipliers (usually in
    (0,1)) to the covariance matrix of the proposal distribution for
    the draws of the intermediate level parameters.  The purpose of this
  \code{Tune} function is to adjust these values so as to achieve
  acceptance ratios of between .2 and .5 in the MCMC draws of the
  \code{THETA}s.}
  
  \item{kappa}{Scalar:  The diagonal of the covariance matrix for the
    (normal) hyperprior distribution for the \eqn{\mu}{mu} parameter.}
  
  \item{nu}{Scalar:  The degrees of freedom for the (Inverse-Wishart)
    hyperprior distriution for the \code{SIGMA} parameter.}
  
  \item{psi}{Scalar:  The diagonal of the matrix parameter of the
    (Inverse-Wishart) hyperprior distribution for the \code{SIGMA}
    parameter.}
  
  \item{mu.vec.0}{Vector:  mean of the (normal) hyperprior distribution for
    the \eqn{\mu}{mu} parameter.}
  
  \item{mu.vec.cu}{Vector of dimension \eqn{R*(C-1)}, where
    \eqn{R}(\eqn{C}) is the number of rows(columns) in each contigency
    table:  Optional starting values for \eqn{\mu}{mu} parameter.}
  
  \item{nolocalmode}{Positive integer:  How often an alternative drawing
  method for the contigency table internal cell counts will be used.
  Use of default value recommended.}
  
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
    contingency table's first row.  Use of default (generated
    internally) recommended.}
  
  \item{numscans}{Positive integer:  How often the algorithm to draw the
  contingency table internal cell counts will be implemented before new
  values of the other parameters are drawn.  Use of default value
  recommended.}
  
  \item{Diri}{Positive integer:  How often a product Dirichlet proposal
    distribution will be used to draw the contingency table row
    probability vectors (the THETAS).}
  
  \item{dof}{Positive integer:  The degrees of freedom of the
    multivariate \eqn{t} proposal distribution used in drawing the contingency
    table row probability vectors (the THETAS).}
  
  \item{debug}{Integer:  Akin to \code{verbose} in some packages.  If set
    to 1, certain status information (including rough notification
    regarding the number of iterations completed) will be
    written to the screen.}
}
\details{
  \code{TuneWithExitPoll} is a necessary precursor function to \code{AnalyzeWithExitPoll}, the workhorse
  function in fitting the R x C
  ecological inference model described in Greiner & Quinn (2009) to
  ecological data augmented by data from surveys of some of the
  contingency tables.  Details and terminology of the basic
  (\emph{i.e.}, without a survey sample) data structure and ecological
  inference model are discussed in the
  documentation accompanying the \code{Analyze} function.  The purpose
  of \code{TuneWithExitPoll}, as prepatory to
  \code{AnalyzeWithExitPoll}, is the same as the purpose of \code{Tune}
  as prepatory to \code{Analyze}.  See the documentation for \code{Tune}
  for a full explanation.

  In the present implementation, the \code{AnalyzeWithExitPoll}, and
  thus \code{TuneWithExitPoll}, presume that the survey consisted of a
  simple random sample from the in-sample contingency tables.
  Future implementations will allow
  incorporation of more complicated survey sampling schemes.

  The arguments to \code{TuneWithExitPoll} are essentially identical
  to those of \code{Tune} with the major exception of
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
  out-of-sample tables), the corresponding row of \code{exitpoll} should
  have a vector of 0s.
}
\value{
  A list with the following elements.
  \item{rhos}{A vector of length \code{I} = number of contingency tables:  each
  element of the \code{rhos} vector is a multiplier used in the proposal
distribution of for draws from the conditional posterior of the THETAs,
as described above.  Feed this vector into the \code{Analyze} function.}

  \item{acc.t}{Matrix of dimension \code{I} x \code{num.runs}:  Each column of
    \code{acc.t} contains the acceptance fractions for the Metropolis-Hastings
    algorithm, with a multivariate \eqn{t_4} proposal distribution, used to draw from the conditional posterior of the
    \code{THETA}s.  If \code{Tune} has worked properly, all elements of the
    final column of this matrix should be between .2 and .5.}

  \item{acc.Diri}{Matrix of dimension \code{I} x \code{num.runs}:  Each column of
    \code{acc.t} contains the acceptance fractions for the Metropolis-Hastings
    algorithm, with independent Dirichlet proposals, used to draw from the conditional posterior of the
    \code{THETA}s.  \code{Tune} does not alter this algorithm.}

  \item{vld.NNs}{A list of length \code{num.runs}:  Each element of
    \code{vld.NNs} is a matrix of dimension \code{I} by \code{R}, with
    each element of the list corresponding to one of the
    \code{num.iters} sets of iterations run by \code{Tune}.  To draw from the conditional posterior of
    the internal cell counts of a contigency table, the \code{Tune} function
    draws R-1 vectors of lenth C from multinomial distributions.  In
    then calculates the counts in the additional row (denote this row as
    r') deterministically.  This procedure can result in negative values
    in row r', in which case the overall proposal for the interior cell
    counts is outside the parameter space (and thus invalid).
    Each matrix of vld.NNs keeps track of the percentage of proposals drawn in
    this manner that are valid (\emph{i.e.}, not invalid).  Each row of
    such a matrix corresponds to a
    contingency table.  Each column in the matrix
    corresponds to a \emph{row} in the a contingency table.  Each entry
    specifies the percentage of multinomial proposals that are valid
    when the specified contingency table row serves as the r' row.  For
    instance, in position 5,2 of vld.NNs is the fraction of valid
    proposals for the 5th contingency table when the second contigency
    table row is the r'th row.  A value of ``NaN'' means that \code{Tune}
    chose to use a different (slower) method of drawing the internal
    cell counts because it suspected that the multinomial method would
    behave badly.}

  \item{acc.NNs}{A list of length \code{num.runs}:  Same as vld.NNs,
    except the entries represent the fraction of proposals accepted (instead of the
    fraction that are in the permissible parameter space).}
}
\references{
  D. James Greiner \& Kevin M. Quinn. 2009.  ``R x C Ecological
  Inference:  Bounds, Correlations, Flexibility, and Transparency of
  Assumptions.'' \emph{J.R. Statist. Soc. A} 172:67-81.}
\author{D. James Greiner, Paul D. Baines, \& Kevin M. Quinn}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
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
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
