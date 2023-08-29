\name{cove.boost.collapse.strata}
\alias{cove.boost.collapse.strata}
\title{
Collapse sample strata for COVE boost correlates study
}
\description{
Collapse sample strata for COVE boost correlates study
}
\usage{

cove.boost.collapse.strata (dat.b, n.demo)

}
\arguments{

  \item{dat.b}{data frame
    }
  \item{n.demo}{
    number of demographics strata, e.g. 6 for COVE correlates
    }
}

\value{
		dat.b, whose Wstratum has been updated
}

\details{
  This function is used by both correlates_processing repo and correlates_reporting3 repo
}
