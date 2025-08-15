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
  \item{data}{A dataframe including individual PFAS measurements.}
  \item{me.PFOSA.AcOH}{The column name of me-PFOSA-AcOH in the dataframe. If there is no measurement of me-PFOSA-AcOH, \code{me.PFOSA.AcOH=NA}. Default is \code{NA}}
  \item{PFUA}{The column name of PFUA in the dataframe. If there is no measurement of PFUA, \code{PFUA=NA}. Default is \code{NA}}
  \item{PFDeA}{The column name of PFDeA in the dataframe. If there is no measurement of PFDeA, \code{PFDeA=NA}. Default is \code{NA}}
  \item{PFNA}{The column name of PFNA in the dataframe. If there is no measurement of PFNA, \code{PFNA=NA}. Default is \code{NA}}
  \item{PFHxS}{The column name of PFHxS in the dataframe. If there is no measurement of PFHxS, \code{PFHxS=NA}. Default is \code{NA}}
  \item{PFOA}{The column name of PFOA in the dataframe. If there is no measurement of PFOA, \code{PFOA=NA}. Default is \code{NA}}
  \item{PFOS}{The column name of PFOS in the dataframe. If there is no measurement of PFOS, \code{PFOS=NA}. Default is \code{NA}}
  
  \item{n.PFOA}{The column name of n-PFOA in the dataframe. If there is no measurement of n-PFOA, \code{n.PFOA=NA}. Default is \code{NA}}
  \item{Sm.PFOS}{The column name of Sm-PFOS in the dataframe. If there is no measurement of Sm-PFOS, \code{Sm.PFOA=NA}. Default is \code{NA}}
  \item{n.PFOS}{The column name of n-PFOS in the dataframe. If there is no measurement of n-PFOS, \code{n.PFOS=NA}. Default is \code{NA}}
  
  \item{use.isomers}{Logical. Whether using the isomers (n-PFOA, Sm-PFOS, and n-PFOS) to calculate PFAS exposure burden or using the sum of isomers (PFOA and PFOS). If \code{use.isomers=TRUE}, the function will calculate PFAS exposure burden using 8 PFAS analytes: me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, n-PFOA, Sm-PFOS, and n-PFOS; if \code{use.isomers=FALSE}, the function will calculate PFAS exposure burden from 7 PFAS analytes: me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, PFOA, and PFOS. Default is \code{FALSE}}
  
  
}
\value{
The function will return a dataframe with individual PFAS measurements, estimated PFAS exposure burden and SE of the estimated PFAS burden.
}
\examples{
set.seed(12345)
data.isomer<-data.frame(matrix(sample(1:100,80,replace=T)/10,ncol=8))

### calculating PFAS exposure burden using isomers
PFAScalc2017(data=data.isomer,
             me.PFOSA.AcOH="X1",
             PFUA="X2",
             PFDeA="X3",
             PFHxS="X4",
             PFNA="X5",
             n.PFOA="X6",
             Sm.PFOS="X7",
             n.PFOS="X8",
             use.isomers=TRUE)
             
### calculating PFAS exposure burden using sum of isomers
data.isomer$X9<-data.isomer$X7+data.isomer$X8

PFAScalc2017(data=data.isomer,
             me.PFOSA.AcOH="X1",
             PFUA="X2",
             PFDeA="X3",
             PFHxS="X4",
             PFNA="X5",
             PFOA="X6",
             PFOS="X9",
             use.isomers=FALSE)

}


