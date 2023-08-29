\name{marginalized.risk.svycoxph.boot}
\alias{marginalized.risk.svycoxph.boot}
\title{Bootstrapping Marginalized Risk Svycoxph Boot}
\description{Functions for plotting.}
\usage{

marginalized.risk.svycoxph.boot(form.0, marker.name, type, data, t, B, 
  ci.type="quantile", numCores=1, additional.terms=NULL, variant=NULL)

}
\arguments{

  \item{form.0}{formula}
  \item{marker.name}{name of the biomarker variable}
  \item{data}{data frame}
  \item{type}{type of risk curves}
  \item{t}{followup time}
  \item{B}{number of bootstrap replicates}
  \item{ci.type}{type of bootstrap confidence interval}
  \item{numCores}{number of parallel processes to run}
  \item{additional.terms}{additional terms to adjust}
  \item{variant}{virus variant}
}

\value{
real
}

\details{
There are three different types.
}
