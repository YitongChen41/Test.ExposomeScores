\name{mod.mirt.sum}
\alias{mod.mirt.sum}
\title{IRT model (summed isomers) for estimating PFAS exposure burden}

\description{
IRT model for estimating PFAS exposure burden calculated from 7 PFAS analytes ( me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, PFOA, and PFOS) in the 2017 - 2018 NHANES. 
}
\examples{
# summary of the IRT model (isomers)
mirt::summary(mod.mirt.sum)

# parameter estimates of the IRT model (summed isomers)
mirt::coef(mod.mirt.sum,simplify=T,IRTpars=T)
}



