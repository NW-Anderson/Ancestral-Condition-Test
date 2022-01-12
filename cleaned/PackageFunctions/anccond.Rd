\name{AncCond}
\alias{AncCond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ancestral Condition Test
}
\description{
Tests whether transitions in a discrete trait are associated with extreme values of a continuous trait.
}
\usage{
AncCondtree, data, drop.state = NULL,
mat = c(0, 2, 1, 0), pi = "estimated",
n.tails = 1, nsim = 100, iter = 100,
message = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
This is the data matrix created by running \code{dataMatrx}
}
## tree: a phylo object
## - for multiphylo analysis. Loop through every tree and compile the the
##   of each loop to create a distribution of observed values and a multiphylo
##   null distribution

## data: a dataframe with three columns (Labels, cont trait, disc trait)
## - Labels should match taxa labels in the phylogeny
## - continuous trait should be numeric values
## - discrete trait must be coded as 1 and 2 if one is ancestral then it must be coded as 1
## There should be no missing data. If a taxa does not have available cont
## and discrete data, please prune it from the tree

## mc: the number of Monte Carlo simulations per simulated dataset
## used to calc p-value

## drop.state: should be NULL unless working under the assumption
## that one state is ancestral and the other derived and back
## transitions are not possible. Using this assumption will
## ignore continuous data from taxa in the derived state

## mat: transition matrix. Should contain the rate
## matrix for evolution of the discrete trait. Acceptable matrices are
## c(0,0,1,0), c(0,1,1,0), c(0,2,1,0)

## pi: The probabilities the root of the tree are either of the
## discrete character states same values possible as make.simmap:
## "equal", "estimated", or vector length 2 with probabilities
## for each state

## n.tails: either 1 or 2 depending on whether user has apriori hypothesis about a certain state

}
\details{
If the rate identity matrix is returned rates are indicated by the numbers:
\cr rate13 demipolyploidy for state2 - odd
\cr
\cr
\cr
The argument \code{constrain} can contain any of the items listed below.  The default state is shown.
\cr
\cr \code{drop.poly=F}
\cr Sets polyploidy rate to zero
\cr

}
\value{
constrained likelihood function is returned
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Heath Blackmon}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}



## INPUT DATA

