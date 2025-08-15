\name{PFAScalc2017}
\alias{hello}
\title{Calculating PFAS exposure burden using the 2017 - 2018 US PFAS exposure burden calculator}
\usage{
PFAScalc2017()
}
\description{
Calculating PFAS exposure burden using the 2017 - 2018 US PFAS exposure burden calculator.
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


