\name{mod.mirt.isomer}
\alias{mod.mirt.isomer}
\title{IRT model (isomers) for estimating PFAS exposure burden}

\description{
IRT model for estimating PFAS exposure burden calculated from 8 PFAS analytes ( me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, n-PFOA, Sm-PFOS, and n-PFOS) in the 2017 - 2018 NHANES. 
}
\examples{
# summary of the IRT model (isomers)
mirt::summary(mod.mirt.isomer)

# parameter estimates of the IRT model (isomers)
mirt::coef(mod.mirt.isomer,simplify=T,IRTpars=T)
}



