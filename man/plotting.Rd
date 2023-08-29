\name{plotting}
\alias{draw.x.axis.cor}
\alias{get.xlim}
\alias{report.assay.values}
\alias{get.marker.histogram}
\title{
Plotting Helper Functions}
\description{
Functions for plotting.
}
\usage{

draw.x.axis.cor(xlim, llox, llox.label, for.ggplot=FALSE)
get.xlim(dat, marker, lloxs) 
report.assay.values(x, assay)
get.marker.histogram(marker, wt, trial, marker.break=marker) 

}
\arguments{

  \item{xlim}{xlim}
  \item{llox}{lower limit}
  \item{llox.label}{label for lower limit}
  \item{for.ggplot}{Boolean}
  \item{dat}{data frame}
  \item{marker}{name of the biomarker variable}
  \item{lloxs}{list of lloxs}
  \item{x}{x}
  \item{assay}{assay}
  \item{wt}{weight}
  \item{trial}{trial}
  \item{marker.break}{marker.break}
}

\value{
real
}

\details{
draw.x.axis.cor is used by both cor_coxph and cor_threshold
}
