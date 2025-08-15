\name{PFAScalc2017}
\alias{PFAScalc2017}
\title{Calculating PFAS exposure burden using the 2017 - 2018 US PFAS exposure burden calculator}
\usage{
PFAScalc2017(data,
  me.PFOSA.AcOH=NA,
  PFUA=NA,
  PFDeA=NA,
  PFHxS=NA,
  PFNA=NA,
  PFOA=NA,
  PFOS=NA,
  n.PFOA=NA,
  Sm.PFOS=NA,
  n.PFOS=NA,
  use.isomers=FALSE)
}
\description{
Calculating PFAS exposure burden using the 2017 - 2018 US PFAS exposure burden calculator.
}
\details{
xxxxxx
}
\arguments{
  \item{data}{A dataframe of individual PFAS measurements.}
  \item{me.PFOSA.AcOH}{The column name of me-PFOSA-AcOH in the dataframe. \code{a}. Can be NA if there is no measurement of me-PFOSA-AcOH. Default is NA. \code{NA}}
}
\value{
The function will return a dataframe with individual PFAS measurements, estimated PFAS exposure burden and SE of the estimated PFAS burden.
}
\examples{
set.seed(12345)
data.isomer<-data.frame(matrix(sample(1:100,80,replace=T)/10,ncol=8))
PFAScalc2017(data=test.isomer,
             me.PFOSA.AcOH="X1",
             PFUA="X2",
             PFDeA="X3",
             PFHxS="X4",
             PFNA="X5",
             n.PFOA="X6",
             Sm.PFOS="X7",
             n.PFOS="X8",
             use.isomers=TRUE)

}


