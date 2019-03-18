\name{Analyze}
\alias{Analyze}

\title{Workhorse Function for Ecological Inference for Sets of R x C Contingency Tables}
\description{
  This function (using the tuned parameters from \code{Tune}) fits a
  hierarchical model to ecological data in which the underlying
  contigency tables can have any number of rows or columns.  The user
  supplies the data and may specify hyperprior values.  Samples from the
  posterior distribution are returned as an \code{mcmc} object, which can be
  analyzed with functions in the \code{coda} package.
}
\usage{
Analyze(fstring, rho.vec, data = NULL, num.iters = 1e+06, 
        save.every =1000, burnin = 10000, 
        mu.vec.0 = rep(log((0.45/(mu.dim - 1))/0.55), mu.dim), 
        kappa = 10, nu = (mu.dim + 6), psi = mu.dim,  
        mu.vec.cu = runif(mu.dim, -3, 0), NNs.start = NULL, 
        THETAS.start = NULL, prob.re = 0.15, sr.probs = NULL, 
        sr.reps = NULL, keep.restart.info = FALSE, 
        keepNNinternals = 0, keepTHETAS = 0, nolocalmode = 50, 
        numscans = 1, Diri = 100, dof = 4, print.every = 10000, 
        debug = 1)
}


\arguments{

  \item{fstring}{String: model formula of contingency tables' column
  totals versus row totals. Must be in specified format (an R
  character string and NOT a true R formula). See Details and
  Examples.}

  \item{rho.vec}{Vector of dimension \eqn{I} = number of contigency
    tables = number of rows in \code{data}: multipliers (usually in
    (0,1)) to the covariance matrix of the proposal distribution for
    the draws of the intermediate level parameters. Typically set to
    the vector output from \code{Tune}.}
 
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

  \item{THETAS.start}{Matrix of dimension \eqn{I} x (\eqn{R*C}),
    where \eqn{I} is the number of contingency tables = number of rows in
    \code{data}:  Optional starting values for the contingency table row
    probability vectors.  The elements in each row of \code{THETAS.start}
    must meet \eqn{R} sum-to-one constraints.  Use of the
    default (randomly generated internally) recommended.}

  \item{prob.re}{A positive fraction:  Probability of random exchange in
  a parallel tempering fitting algorithm.  Not yet implemented.}

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

  \item{print.every}{Positive integer:  If \code{debug == 1}, the number
    of every \code{print.every}\eqn{th} iteration will be written to the
    screen.  Must be divisible into \code{num.iters}.}

  \item{debug}{Integer:  Akin to \code{verbose} in some packages.  If set
    to 1, certain status information (including rough notification
    regarding the number of iterations completed) will be
  written to the screen.}

}

\details{
  \code{Analyze} is the workhorse function in fitting the R x C
  ecological inference model described in Greiner & Quinn (2009).
  
  Ecological data consist of sets of contingency tables in which the row
  and column totals, but none of the internal cell counts, are
  observed.  For example, in the context of voting rights litigation,
  there is often one contigency table for each voting precinct; the row
  totals are voting-age population figures, with each row representing a
  race/ethnicity; all but the last (right-most) column representing
  votes cast for particular candidates; and the last (right-most) column
  representing persons not voting.

  The model described in Greiner & Quinn (2009) conditions on the row
  totals throughout. In each contigency table, the rows are assumed to
  follow mutually independent multinomials, conditional on separate
  probability vectors which are denoted \eqn{theta_r}{\theta_r} for
  \eqn{r = 1} to \eqn{R} (\eqn{R} being the number of rows in each
  contigency table). Each \eqn{theta_r}{\theta_r} then undergoes a
  multidimensional logistic transformation, using the last
  (right-most) column as the reference category. This results in \eqn{R}
  transformed 
  vectors of dimension \eqn{(C-1)}; these transformed vectors, denoted
  \eqn{\omega_rs}{omega_r's}, are stacked to form a single
  \eqn{\omega}{omega} vector corresponding to that contingency
  table. The omega vectors are assumed to follow (i.i.d.)  a
  multivariate normal distribution. A standard N(\eqn{\mu_0}{mu},
  \eqn{\kappa}{kappa} * I) and Inv-Wish(\eqn{\nu}{nu}, \eqn{\psi}{psi}
  * I) (I is the identity matrix) prior is placed on the normal. The
  user may set \eqn{\mu_0}{mu}, \eqn{\kappa}{kappa}, \eqn{\nu}{nu},
  amd \eqn{\psi}{psi}.

  fstring must be in a specific format.  It must be a string, and it
  must consist of (i) the names of vectors of contingency table column
  totals separated by commas, (ii) then a tilde, (iii) then the names of vectors of
  contingency table row
  totals separated by commas. The order in which the contigency table
  column totals are listed is important because the final column
  with become the reference category in the multidimensional logistic
  transformation described above.  See Examples.

  Fitting the model is accomplished via a Gibbs sampler in which the
  internal cell counts (for each contingency table), then the
  \eqn{\theta's}{thetas}, and then the \eqn{\mu}{mu} and
  \eqn{\Sigma}{Sigma} parameters are drawn in turn.  This method
  automatically produces draws of the internal cell counts, functions of
  which are often the true targets of inference.

  The function returns an object of class mcmc suitable for use in
  functions from the coda package, including combination
  (with other outputs from Analyze) into an object of class mcmc.list.
  The return object includes draws from the posterior distribution of
  the following items:  each element of the \eqn{\mu}{mu}; the standard
  deviations in \eqn{\Sigma}{Sigma} (meaning the square root of the
  diagonal elements); the correlations in \eqn{\Sigma}{Sigma}; the
  sums across all contigency tables of the counts in each of the R * C
  internal cell positions; and a series of functions of these sums that
  are often of interest in voting applications (these may obviously be
  ignored if interest lies elsewhere).  Except for the correlations from
  \eqn{\Sigma}{Sigma}, the labeling follows a self-evident pattern, with
  the names taking from \code{fstring}.  The correlations are
  labeled by a combination of two numbers, representing their position
  in the \eqn{\Sigma}{Sigma} matrix.

  The series of functions of the internal cell counts calculated
  automatically fall into four categories: \code{LAMBDA},
  \code{TURNOUT}, \code{GAMMA}, and \code{Beta}. To explain
  these terms, consider an example in which the contigency tables have
  three rows ("bla", "whi", and "his") and three columns ("Dem",
  "Rep", "Abs"), corresponding to black, white, Hispanic, Democratic,
  Republican, and Abstain (from voting). Thus, in position 1,1 of each
  contigency table is the (unobserved) number of blacks voting
  Democrat, position 2,1 holds the (unobserved) number of whites
  voting Democrat, etc. In each position (1,1 = black Democrat votes;
  2,1 = white Democrat votes), sum across all \eqn{I} contingency
  tables to produce a single R x C table consisting of summed counts
  (these sums are, incidentally, the NN values the software reports).
  \code{LAMBDA}, \code{TURNOUT}, \code{GAMMA}, and
  \code{Beta} are functions of these summed counts, as explained
  in the paragraph below. Note that the paragraph below refers to the
  counts in the single table produced by this summing process.
  Notation: \eqn{NN_{rc}}{NN_rc} is the sum (over all contingency
  tables) of the counts in cell r,c. So \eqn{NN_{bD}}{NN_bD} is the
  total number of blacks voting Democract, \eqn{NN_{wD}}{NN_wD} is the
  total number of whites voting Democrat, etc.

  \code{LAMBDA}: For example, \code{LAMBDA_bD} = \eqn{NN_{bD}/(NN_b -
  NN_{bA})}{NN_bD/(NN_b - NN_bA)}.  Similarly, \code{LAMBDA_hR} = 
  \eqn{NN_{hR}/(NN_{h} - NN_{hA})}{NN_hR/(NN_h - NN_hA)}. In voting
  parlance, this the fraction of each race's \emph{voters} supporting
  a particular candidate. There are \eqn{R * (C-1)} LAMBDAs, \eqn{C-1}
  of them for each row.

  \code{TURNOUT}: For example, \eqn{Turnout_w = (NN_{w} -
    NN_{wA})/NN_w}{Turnout_w = (NN_w - NN_wA)/NN_w}. In voting 
  parlance, this is the fraction of each race that showed up to vote.
  There are R TURNOUTs, one for each row.

  \code{GAMMA}: For example, \code{GAMMA_h} = \eqn{(NN_{hD} +
  NN_{hR})/(NN_{bD} + NN_{bR} + NN_{wD} + NN_{wR} + NN_{hD} +
  NN_{hR})}{(NN_hD + NN_hR)/(NN_bD + NN_bR + NN_wD + NN_wR + NN_hD +
  NN_hR)}. In voting parlance, this is the fraction that each race
  contributes to the voting electorate.

 \code{BETA}: For example, \code{BETA_wR} = \eqn{(NN_{wR})/(NN_{wD} +
  NN_{wR} + NN_{wA}) = (NN_{wR})/(NN_w)}{(NN_wR)/(NN_wD + NN_wR + NN_wA) =
  (NN_wR)/(NN_w)}. This is the fraction of each race's potential (as
  opposed to actual) voters supporting a
  particular candidate. Although there are theoretically \eqn{R * C}
  \code{BETA} values that could be calculated, in fact the
  \code{BETA} values for the last (reference) category are ignored,
  so only \eqn{R * (C-1)} are calculated.

  If \code{keepNNinternals} is non-zero, the specified number of draws
  of the internal cell counts for each contingency table will be save.
  These may be retreived via \code{attr} (see Examples, below).  The
  result is a matrix of dimension \code{keepNNinternals} \eqn{\times R
  * C * I}{x R * C * I}, where \eqn{I} is the number of contigency
  tables.  Each row consists of an iteration's draws.  The first
  column contains the draws of the counts in position 1,1 in the first
  contigency table, the second colum contains the draws of position
  1,2 in the first contigency table, etc.  In other words, the columns
  in the output represent the first contigency table vectorized row
  major, then the second contingency table vectorized row major, etc.
  The same applies to \code{keepTHETAS}, except that the \code{THETA}s
  represent the multinomial row probabilties.}

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

\section{Warnings}{
  
  \strong{Computer time:} At present, using this function (and the others
    in this package) requires substantial computer time.  The lack of
    information in ecological data results
    in slow mixing chains, and the number of parameters that must be drawn
    in each Gibbs sampler iteration is large.  Chain length should be
    adjusted to achieve adequate convergence.  In general, the more
    segregated the housing patterns in the jurisdiction (meaning the
    greater the percentage of contingency tables in which one row's counts
    make up a large portion of that table's total), the smaller the number of
    iterations needed.  We are exploring more efficient
    sampling algorithms that we anticipate will result in better mixing
    and faster drawing.  At present, however, users should anticipate that
    analysis of a dataset will take several hours.
 
  \strong{Large datasets:} At present, use of this fuction (and thus
    this package) is not recommended for large (\emph{i.e.}, more than
    1000 contingency tables) datasets. See immediately above.
     
  \strong{RAM requirements:} Do not select large values of
    \code{keepNNinternals} or \code{keepTHETAS} without adequate RAM.

  \strong{Gelman-Rubin diagnostic in the CODA package:} Using the
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
library(RxCEcolInf)
data(stlouis)
Tune.stlouis <- Tune("Bosley, Roberts, Ribaudo, Villa, NoVote ~ bvap, ovap",
                     data = stlouis,
                     num.iters = 10000,
                     num.runs = 15)
Chain1.stlouis <- Analyze("Bosley, Roberts, Ribaudo, Villa, 
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
names(attributes(stlouis.MCMClist))
summary(stlouis.MCMClist, quantiles = c(.025, .05, .5, .95, .975))
plot(stlouis.MCMClist)
geweke.diag(stlouis.MCMClist)
heidel.diag(stlouis.MCMClist)
#  Do not run gelman.diag; see warnings
NNs <- attr(stlouis.MCMClist, "NN.internals")
THETAS <- attr(stlouis.MCMClist, "THETA")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% __ONLY ONE__ keyword per line

