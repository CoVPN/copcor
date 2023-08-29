\name{utils}
\alias{marker.name.to.assay}
\alias{get.marginalized.risk.no.marker}
\alias{add.trichotomized.markers}
\alias{make.case.count.marker.availability.table}
\title{
Utility Functions
}
\description{
Helpful functions.
}
\usage{

marker.name.to.assay (a)
get.marginalized.risk.no.marker(formula, dat.ph1, followup.day)
add.trichotomized.markers(dat, markers, ph2.col.name="ph2", wt.col.name="wt") 
make.case.count.marker.availability.table(dat)

}
\arguments{

  \item{a}{assay name}
  \item{formula}{formula}
  \item{dat.ph1}{phase 1 data frame}
  \item{dat}{data frame}
  \item{followup.day}{follow up day for computing risk}
  \item{markers}{strings}
  \item{ph2.col.name}{name of the ph2 indicator column}
  \item{wt.col.name}{name of the weight column}
}

\value{
		string
}

\details{
  This function ...
}
