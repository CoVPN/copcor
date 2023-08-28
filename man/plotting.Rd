\name{plotting}
\alias{draw.x.axis.cor}
\alias{get.xlim}
\title{
Plotting Helper Functions}
\description{
Functions for plotting.
}
\usage{

draw.x.axis.cor(xlim, llox, llox.label, for.ggplot=FALSE)
get.xlim(dat, marker, lloxs) 

}
\arguments{

  \item{xlim}{xlim
    }
  \item{llox}{lower limit
    }
  \item{llox.label}{label for lower limit
    }
  \item{for.ggplot}{Boolean    }
  \item{dat}{data frame}
  \item{marker}{name of the biomarker variable}
  \item{lloxs}{list of lloxs}
}

\value{
real
}

\details{
draw.x.axis.cor is used by both cor_coxph and cor_threshold
}
