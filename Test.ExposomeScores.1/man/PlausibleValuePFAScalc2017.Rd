\name{PlausibleValuePFAScalc2017}
\alias{PlausibleValuePFAScalc2017}
\title{Calculating plausible values of PFAS exposure burden}
\usage{
PlausibleValuePFAScalc2017(data,
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
  ndraws=50,
  seed=12345,
  use.isomers=FALSE
  )
}
\description{
Calculating plausible values of PFAS exposure burden using the 2017 - 2018 US PFAS exposure burden calculator. These plausible values may enable users to verify robustness of PFAS-outcome associations while accounting for standard errors of PFAS burden scores.
}

\arguments{
  \item{data}{A dataframe including individual PFAS measurements. All PFAS concentrations in the data need to be measured in ng/mL.}
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
  
  \item{ndraws}{Number of plausible values for each individual. Default is \code{ndraws=5}}
  \item{seed}{Specify random seed. Default is \code{seed=12345}}
  
  \item{use.isomers}{Logical. Whether using the isomers (n-PFOA, Sm-PFOS, and n-PFOS) to calculate PFAS exposure burden or using the sum of isomers (PFOA and PFOS). If \code{use.isomers=TRUE}, the function will calculate PFAS exposure burden using 8 PFAS analytes: me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, n-PFOA, Sm-PFOS, and n-PFOS; if \code{use.isomers=FALSE}, the function will calculate PFAS exposure burden from 7 PFAS analytes: me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, PFOA, and PFOS. Default is \code{FALSE}}
  
  
}

\details{
This function calculates the plausible values of estimated PFAS burden scores basing on two IRT models which use survey-weighted decile cutoffs for PFAS analytes calculated from biomonitoring data from 1915 children and adults in the 2017-2018 NHANES. (Liu et al., 2022)

PFAS burden (using isomers of PFOA and PFOS) is calculated from 8 PFAS analytes: me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, n-PFOA, n-PFOS, and Sm-PFOS. 

PFAS burden (using summed PFOA and PFOS) is calculated from 7 PFAS analytes: me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, PFOA, and PFOS. In this model, we use the sum concentrations of isomers for PFOA and PFOS instead of using the concentrations of branched and linear isomers (n-PFOS, Sb-PFOA, n-PFOS, and Sm-PFOS).

\strong{Plausible values of the estimated PFAS exposure burden score)} \cr
The PFAS exposure burden score estimated from IRT models have standard errors. This function can calculate plausible values of the estimated PFAS burden for researchers to perform secondary analyses to account for the uncertainty in latent burden scores. 
}

\references{
Liu, S. H., J. R. Kuiper, Y. Chen, L. Feuerstahler, J. Teresi, and J. P. Buckley. 2022. 'Developing an Exposure Burden Score for Chemical Mixtures Using Item Response Theory, with Applications to PFAS Mixtures', Environ Health Perspect, 130: 117001. \url{https://ehp.niehs.nih.gov/doi/10.1289/EHP10125}
}


\value{
The function will return a list of plausible values for the estimated PFAS exposure burden.
}
\examples{
set.seed(12345)
plausible.values<-PlausibleValuePFAScalc2017(data=PFAScycleJ,
                          me.PFOSA.AcOH="me.PFOSA.AcOH",
                          PFUA="PFUA",
                          PFDeA="PFDeA",
                          PFHxS="PFHxS",
                          PFNA="PFNA",
                          n.PFOA="n.PFOA",
                          Sm.PFOS="Sm.PFOS",
                          n.PFOS="n.PFOS",
                          ndraws=5,
                          seed=12345,
                          use.isomers=TRUE)

# 5 plausible values of estimated PFAS burden scores for each participant
head(do.call("cbind",plausible.values))


}


