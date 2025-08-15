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

\details{
This function calculates the PFAS burden scores basing on two IRT models which use survey-weighted decile cutoffs for PFAS biomarkers calculated from biomonitoring data from 1915 children and adults in the 2017-2018 NHANES. (Liu et al., 2022)

PFAS burden (using isomers of PFOA and PFOS) is calculated from 8 PFAS biomarkers: me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, n-PFOA, n-PFOS, and Sm-PFOS. 

PFAS burden (using summed PFOA and PFOS) is calculated from 7 PFAS biomarkers: me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, PFOA, and PFOS. In this model, we use the sum concentrations of isomers for PFOA and PFOS instead of using the concentrations of branched and linear isomers (n-PFOS, Sb-PFOA, n-PFOS, and Sm-PFOS).

All PFAS biomarker concentrations in the data need to be measured in ng/mL

\strong{Quantifying PFAS exposure burden scores using item response theory (IRT)} \cr

We used IRT to understand how the set of assayed PFAS biomarkers provides information about the underlying latent exposure burden to PFAS. IRT-based scores represent estimates of the underlying (latent) PFAS burden and differentially weight the contributions of the individual biomarkers based on their data-driven difficulty and discrimination parameters, with a nonlinear relationship between levels of individual biomarkers and the latent burden.

IRT is a flexible modeling approach with several similarities to confirmatory factor analysis (CFA), which attempts to describe and elucidate the relations of measured items (e.g. PFAS biomarker levels) as indicators for an unobservable latent construct (e.g. latent exposure burden).

In an IRT framework, the probability of an individual endorsing an item (e.g. higher PFOA level) is dependent on several factors, including the latent construct (e.g. PFAS burden) and item-specific parameters. Additionally, the probability of endorsing an item increases monotonically with the level of the trait.


\strong{For more information, please refer to our paper.} \cr

https://ehp.niehs.nih.gov/doi/10.1289/EHP10125

Liu, S. H., J. R. Kuiper, Y. Chen, L. Feuerstahler, J. Teresi, and J. P. Buckley. 2022. 'Developing an Exposure Burden Score for Chemical Mixtures Using Item Response Theory, with Applications to PFAS Mixtures', Environ Health Perspect, 130: 117001.

}

\references{https://ehp.niehs.nih.gov/doi/10.1289/EHP10125}


\value{
The function will return a dataframe with individual PFAS measurements, estimated PFAS exposure burden and SE of the estimated PFAS burden.
}
\examples{
set.seed(12345)
data.isomer<-data.frame(matrix(sample(1:100,80,replace=T)/10,ncol=8))

# calculating PFAS exposure burden using isomers
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

# calculating PFAS exposure burden using sum of isomers        
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


