\name{gendata.ep}
\alias{gendata.ep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function To Simulate Ecological and Survey Data For Use in
  Testing And Analyzing Other Functions in Package}
\description{
  This function generates simulated ecological data, \emph{i.e.}, data in the form of contigency
  tables in which the row and column totals but none of the internal
  cell counts are observed.  At the user's option, data from simulated surveys of
  some of the `units' (in voting parlance, 'precincts') that gave rise
  to the contingency tables are also produced.
}

\usage{ 
gendata.ep(nprecincts = 175,
           nrowcat = 3,
           ncolcat = 3,
           colcatnames = c("Dem", "Rep", "Abs"),
           mu0 = c(-.6, -2.05, -1.7, -.2, -1.45, -1.45),
           rowcatnames = c("bla", "whi", "his", "asi"),
           alpha = c(.35, .45, .2, .1),
           housing.seg = 1,
           nprecincts.ep = 40,
           samplefrac.ep = 1/14,
           K0 = NULL,
           nu0 = 12,
           Psi0 = NULL,
           lambda = 1000,
           dispersion.low.lim = 1,
           dispersion.up.lim = 1,
           outfile=NULL,
           his.agg.bias.vec = c(0,0),
           HerfInvexp = 3.5,
           HerfNoInvexp = 3.5,
           HerfReasexp = 2)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nprecincts}{positive integer:  The number of contingency tables (precincts)
    in the simulated dataset.}
  
  \item{nrowcat}{integer > 1:  The number of rows in each of the
    contingency tables.}
  
  \item{ncolcat}{integer > 1:  The number of columns in each of the
    contingency tables.}

  \item{rowcatnames}{string of length = length(\code{nrowcat}):  Names
    of rows in each contingency table.}
  
  \item{colcatnames}{string of length = length(\code{ncolcat}):  Names
    of columns in each contingency table.}
  
  \item{alpha}{vector of length(\code{nrowcat}):  initial parameters to
    a  Dirichlet distribution used to generate each contingency table's
    row fractions.}
  
  \item{housing.seg}{scalar > 0:  multiplied to alpha to generate final
    parameters to Dirichlet distribution used to generate each
    contingency table's row fractions.}
  
  \item{mu0}{vector of length (\code{nrowcat} * (\code{ncolcat} - 1)):
    The  mean of the
    multivariate normal hyperprior at the top level of the hierarchical
    model from which the data are simulated.  See Details.}
  
  \item{K0}{square matrix of dimension (\code{nrowcat} * (\code{ncolcat}
    - 1)):  the covariance
    matrix of the multivariate normal hyperprior at the top level of the
  hierarchical model from which the data are simulated.  See Details.}
  
  \item{nu0}{scalar > 0:  the degrees of freedom for the Inv-Wishart
    hyperprior from which the \eqn{SIGMA}{\Sigma} matrix will be drawn.}
  
  \item{Psi0}{square matrix of dimension (\code{nrowcat} *
    (\code{ncolcat} - 1)):  scale matrix for the Inv-Wishart
    hyperprior from which the \code{SIGMA} matrix will be drawn.}
  
  \item{lambda}{scalar > 0:  initial parameter of the Poisson
    distribution from which the number of voters in each precinct will be drawn}
  
  \item{dispersion.low.lim}{scalar > 0 but < dispersion.up.lim:
    lower limit of a draw from \code{runif()} to be multiplied to
    \code{lambda} to set a lower limit on
    the parameter used to draw from the Poisson distribution that
    determines the number of voters in each precinct.}
  
  \item{dispersion.up.lim}{scalar > dispersion.low.lim:
    upper limit of a draw from \code{runif()} to be multiplied
    to \code{lambda} to set a upper limit on
    the parameter used to draw from the Poisson distribution that
    determines the number of voters in each precinct.}
  
  \item{outfile}{string ending in ".Rdata":  filepath and name of
    object; if non-NULL, the object returned by this function will be
    saved to the location specified by \code{outfile}.}
  
  \item{his.agg.bias.vec}{vector of length 2:  only implemented for
    nowcat = 3 and ncolcat = 3:  if non-null, induces aggregation bias
    into the simulated data.  See Details.}
  
  \item{nprecincts.ep}{integer > -1 and less than nprecincts:  number of
    contingency tables (precincts) to be included in simulated survey
    sample (ep for "exit poll").}
  
  \item{samplefrac.ep}{fraction (real number between 0 and 1):
    percentage of individual units (voters) within each contingency table
    (precinct) include in the survey sample.}
  
  \item{HerfInvexp}{scalar:  exponent used to generate inverted
    quasi-Herfindahl weights used to sample contingency tables (precincts) for
  inclusion in a sample survey.  See Details.}
  
  \item{HerfNoInvexp}{scalar:  same as HerInvexp except the quasi-Herfindahl
    weights are not inverted.  See Details.}
  
  \item{HerfReasexp}{scalar:  same as HerfInvexp, for a separate sample
    survey.  See Details.}
}
\details{
  This function simulates data from the ecological inference model outlined in Greiner \&
  Quinn (2009).  At the user's option (by setting nprecincts.ep to an
  integer greater than 0), the function generates three survey samples
  from the simulated dataset.  The specifics of the function's operation
  are as follows.

  First, the function simulates the total number of individual units
  (voters) in each contigency table (precinct) from a Poisson
  distribution with parameter \code{lambda} * runif(1, dispersion.low.lim,
  dispersion.up.lim).  Next, for each table, the function simulates the
  vector of fraction of units (voters) in each table (precinct) row.
  The fractions are simulated from a Dirichlet distribution with
  parameter vector \code{housing.seg} * \code{alpha}.  The row fractions are
  multiplied by the total number of units (voters), and the resulting
  vector is rounded to produce contingency table row counts for each
  table.

  Next, a vector \eqn{mu}{\mu} is simulated from a multivariate normal
  with mean \code{mu0} and covariance matrix \code{K0}.  A covariance
  matrix \code{Sigma} is simulated from an Inv-Wishart with
  \code{nu0} degrees of freedom and scale matrix \code{Psi0}.

  Next, \code{nprecincts} vectors are drawn from \eqn{N(mu,
    SIGMA)}{N(\mu, \Sigma)}.  Each of these draws undergoes an inverse-stacked
  multidimensional logistic transformation to produce a set of \code{nrowcat}
  probability vectors (each of which sums to one) for \code{nrowcat}
  multinomial distributions, one for each row in that contingency
  table.  Next, the \code{nrowcat} multinomial values, which represent the true (and
  in real life, unobserved) internal cell counts, are drawn from the relevant row
  counts and these probability vectors.  The column totals are
  calculated via summation.

  If \code{nprecincts.ep} is greater than 0, three simulated surveys (exit polls) are
  drawn.  All three select contingency tables (precincts) using weights
  that are a function of the composition of the row totals.  Specifically the row
  fractions are raised to a power q and then summed (when q = 2 this calculation is
  known in antitrust law as a Herfindahl index).  For one of the three
  surveys (exit polls) \code{gendata.ep} generates, these
  quasi-Herfindahl indices are the weights.  For two of the three
  surveys (exit polls) \code{gendata.ep} generates, denoted \code{EPInv}
  and \code{EPReas}, the sample weights are the reciprocals of these
  quasi-Herfindhal indices.  The former method tends to weight
  contingency tables (precincts) in which one row dominates the table
  higher than contigency tables (precincts) in which row fractions are close to the
  same.  In voting parlance, precincts in which one racial group
  dominates are more likely to be sampled than racially mixed
  precincts.  The latter method, in which the sample weights are
  reciprocated, weights contingency tables in which row fractions are
  similar more highly; in voting parlance, mixed-race precincts are more
  likly to be sampled.

  For example, suppose \code{nrowcat} = 3, \code{HerInvexp} = 3.5,
  \code{HerfReas} = 2, and
  \code{HerfNoInv} = 3.5.  Consider
  contingency table P1 with row counts (300, 300, 300) and contingency
  table P2 with row counts (950, 25, 25).  Then:

  \bold{Row fractions:} The corresponding row
  fractions are (300/900, 300/900, 300/900) = (.33, .33, .33) and
  (950/1000, 25/1000, 25/1000) = (.95, .025, .025).

  \bold{EPInv weights:} \code{EPInv} would
  sample from assign P1 and P2 weights as follows:  \eqn{1/sum(.33^3.5,
  .33^3.5, .33^3.5) = 16.1} and \eqn{1/sum(.95^3.5, .025^3.5, .025^3.5) =
  1.2}.

  \bold{EPReas weights:} \code{EPReas} would assign weights as
  follows: \eqn{1/sum(.33^2, .33^2, .33^2) = 3.1} and \eqn{1/sum(.95^2, .025^2,
  .025^2) = 1.1}.

  \bold{EPNoInv weights:} \code{EPNoInv} would assign weights as
  follows: \eqn{sum(.33^3.5, .33^3.5, .33^3.5) = .062} and \eqn{sum(.95^3.5,
  .025^3.5, .025^3.5) = .84}.

  For each of the three simulated surveys (\code{EPInv}, \code{EPReas},
  and \code{EPNoInv}), \code{gendata.ep} returns a list of length
  three.  The first element of the list, \code{returnmat.ep}, is a matrix of
  dimension \code{nprecincts} by (\code{nrowcat} * \code{ncolcat})
  suitable for passing to \code{TuneWithExitPoll} and
  \code{AnalyzeWithExitPoll}.  That is, the first row of
  \code{returnmat.ep} corresponds to the first row of \code{GQdata},
  meaning that they both contain information from the same
  contingency table.  The second row of \code{returnmat.ep} contains
  information from the contingency table represented in the second row
  of \code{GQdata}.  And so on.  In addition, \code{returnmat.ep} has counts
  from the sample of the contingency table in vectorized row major
  format, as required for \code{TuneWithExitPoll} and
  \code{AnalyzeWithExitPoll}.

  If \code{nrowcat} = \code{ncolcat} = 3, then the user may set
  \code{his.agg.bias.vec} to be nonzero.  This will introduce aggregation
  bias into the data by making the probability vector of the second row
  of each contingency table a function of the fractional composition of
  the third row.  In voting parlance, if the rows are black, white, and
  Hispanic, the white voting behavior will be a function of the percent
  Hispanic in each precinct.  For example, if \code{his.agg.bias.vec} =
  c(1.7, -3), and if the fraction Hispanic in each precinct i is
  \eqn{X_hi}{X_{h_i}}, then in the ith precinct, the \eqn{mu_i[3]}{\mu_i[3]}
  is set to \code{mu0[3]} + \eqn{X_hi * 1.7}{X_{h_i} * 1.7}, while \eqn{mu_i[4]}{\mu_i[4]}
  is set to \code{mu0[4]} + \eqn{X_hi * -3}{X_{h_i} * -3}.  This feature
  allows testing of the ecological inference model with aggregation
  bias.
  
}
\value{
  A list with the follwing elements.
  \item{GQdata}{Matrix of dimension \code{nprecincts} by (\code{nrowcat}
  + \code{ncolcat}):  The simulated (observed) ecological data, meaning
  the row and column totals in the contingency tables.  May be passed as
  \code{data} argument in \code{Tune}, \code{Analyze},
  \code{TuneWithExitPoll}, and \code{AnalyzeWithExitPoll}}
  
  \item{EPInv}{List of length 3:  \code{returnmat.ep}, the
  first element in the list, is a matrix that may be passed as the
  \code{exitpoll} argument in \code{TuneWithExitPoll} and
  \code{AnalyzeWithExitPoll}. See Details.  \code{ObsData} is
  a dataframe that may be used as the \code{data} argument in the
  \code{survey} package.  \code{sampprecincts.ep} is a vector detailing
  the row numbers of \code{GQdata} (meaning the contingency tables) that
  were  included in the \code{EPInv} survey (exit
  poll).  See Details for an explanation of the weights used to select the contingency
  tables for inclusion in the \code{EPInv} survey (exit poll).}

  \item{EPNoInv}{List of length 3:  Contains the same elements as
    \code{EPInv}.  See Details for an explanation of weights used to
    select the contingency tables for inclusion in the \code{EPNoInv}
    survey (exit poll).}

  \item{EPReas}{List of length 3:  Contains the same elements as
    \code{EPInv}.  See Details for an explanation of weights used to
    select the contingency tables for inclusion in the \code{EPReas}
    survey (exit poll).}

  \item{omega.matrix}{Matrix of dimension \code{nprecincts} by
    (\code{nrowcat} * (\code{ncolcat}-1)):  The matrix of draws from the
    multivariate normal distribution at the second level of the hiearchical
    model giving rise to \code{GQdata}.  These values undergo an
    inverse-stacked-multidimensional logistic transformation to produce contingency
    table row probability vectors.}

  \item{interior.tables}{List of length \code{nprecincts}:  Each element of the
    list is a full (meaning all interior cells are filled in)
    contingency table.}

  \item{mu}{vector of length \code{nrowcat} * (\code{ncolcat}-1):  the
    \eqn{mu}{\mu} vector drawn at the top level of the hierarchical
    model giving rise to \code{GQdata}.  See Details.}

  \item{Sigma}{square matrix of dimension \code{nrowcat} *
    (\code{ncolcat}-1): the covariance matrix  drawn at the top level of
    the hierarchical model giving rise to \code{GQdata}.  See Details.}

  \item{Sigma.diag}{the output of \code{diag(Sigma)}.}

  \item{Sigma.corr}{the output of \code{cov2cor(Sigma)}.}

  \item{sim.check.vec}{vector:  the true values of the
    parameters generated by \code{Analyze} and
    \code{AnalyzeWithExitPoll} in the same order as the parameters are
    produced by those two functions.  This vector is useful in assessing
    the coverage of intervals from the posterior draws from \code{Analyze} and
    \code{AnalyzeWithExitPoll}.}
}
\references{D. James Greiner \& Kevin M. Quinn. 2009.  ``R x C Ecological
  Inference:  Bounds, Correlations, Flexibility, and Transparency of
  Assumptions.'' \emph{J.R. Statist. Soc. A} 172:67-81.}
\author{D. James Greiner \& Kevin M. Quinn}

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