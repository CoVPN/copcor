\name{bootstrap_utils}
\alias{get.bootstrap.data.cor}
\alias{bootstrap.cove.boost}
\alias{bootstrap.cove.boost.2}
\alias{bootstrap.case.control.samples}
\alias{get.ptids.by.stratum.for.bootstrap}
\title{
Bootstrap Utility Functions
}
\description{
Helpful functions.
}
\usage{

get.bootstrap.data.cor (data, ptids.by.stratum, seed) 
bootstrap.cove.boost (dat.ph1, seed) 
bootstrap.cove.boost.2 (dat.ph1, seed) 
bootstrap.case.control.samples (dat.ph1, seed, delta.name="EventIndPrimary", 
  strata.name="tps.stratum", ph2.name="ph2", min.cell.size=1) 
get.ptids.by.stratum.for.bootstrap (data) 

}
\arguments{
  \item{data}{data frame}
  \item{seed}{rng seed}
  \item{dat.ph1}{phase 1 data frame}
  \item{ptids.by.stratum}{a list}
  
  \item{delta.name}{name of event indicator column}
  \item{strata.name}{name of the strata indicator column}
  \item{ph2.name}{name of the ph2 indicator column}
  \item{min.cell.size}{minimum required cell size}
}

\value{
		string
}

\details{
  This function ...
}
