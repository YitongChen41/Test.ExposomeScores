\name{PFASburdenNHANEScomparison}
\alias{PFASburdenNHANEScomparison}
\title{Plotting the distribution of PFAS exposure burden scores between sample population and the age- and sex-matched participants in NHANES}
\usage{
PFASburdenNHANEScomparison(ID,Sex,Age,PFAS.burden,use.isomers=FALSE,matched.ratio=5,bins.num=50,plot.density=TRUE)
}
\description{
This function plots the distribution (density or histogram) of the PFAS exposure burden scores between the study sample population and the age- and sex-matched participants in the 2017-2018 NHANES cycle
}

\arguments{
  \item{ID}{The unique identifier of the study sample.}
  \item{Sex}{A factor of participants' sex. 1=Male; 0=Female.}
  \item{Age}{A numeric vector of participants' age in years.}
  \item{PFAS.burden}{Estimated PFAS exposure burden of the study sample.}
  \item{use.isomers}{Logical. Whether using the isomers (n-PFOA, Sm-PFOS, and n-PFOS) to calculate PFAS exposure burden or using the sum of isomers (PFOA and PFOS). If \code{use.isomers=TRUE}, the estimated PFAS exposure burden of the study sample is compared with the PFAS exposure burden of the matched NHANES sample calculated from 8 PFAS analytes (me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, n-PFOA, Sm-PFOS, and n-PFOS); if \code{use.isomers=FALSE}, the estimated PFAS exposure burden of the study sample is compared with the PFAS exposure burden of the matched NHANES sample calculated from 7 PFAS analytes (me-PFOSA-AcOH, PFUA, PFDeA, PFHxS, PFNA, PFOA, and PFOS). Default is \code{FALSE}}
  \item{matched.ratio}{The matched.ratio indicates how many participants from NHANES should be matched (with replacement) to each individual in the study sample. It should be a single integer value. The default is \code{matched.ratio=5}. Each individual in the study sample will be matched with 5 participants in the 2017-2018 NHANES cycle, with exact matching by sex and neareast matching by age.}
  \item{bins.num}{Number of bins in generating the histogram. Default is \code{bins.num=50}}
  \item{plot.density}{Logical. If \code{plot.density=TRUE}, the function will draw the kernel density estimates, which is a smoothed version of the histogram; if \code{plot.density=FALSE}, the function will draw the histogram. Default is \code{plot.density=TRUE}}
}

\examples{

# participants from the 2015-2016 NHANES cycle
head(PFAScycleI)

# calculating PFAS exposure burden using sum of isomers        
EstimatedPFASburden<-PFAScalc2017(data=PFAScycleI,
             me.PFOSA.AcOH="me.PFOSA.AcOH",
             PFUA="PFUA",
             PFDeA="PFDeA",
             PFHxS="PFHxS",
             PFNA="PFNA",
             PFOA="PFOA",
             PFOS="PFOS",
             use.isomers=FALSE)$pfas.burden$`PFAS burden`

# Histogram of estimated PFAS burden score compared with the sex- and age-matched participants in 2017-2018 NHANES cycle
PFASburdenNHANEScomparison(ID=PFAScycleI$ID,
                           Sex=PFAScycleI$Gender,
                           Age=PFAScycleI$Age,
                           PFAS.burden=EstimatedPFASburden,
                           use.isomers=FALSE,
                           matched.ratio=1,
                           plot.density=FALSE)

# Density of estimated PFAS burden score compared with the sex- and age-matched participants in 2017-2018 NHANES cycle
PFASburdenNHANEScomparison(ID=PFAScycleI$ID,
                           Sex=PFAScycleI$Gender,
                           Age=PFAScycleI$Age,
                           PFAS.burden=EstimatedPFASburden,
                           use.isomers=FALSE,
                           matched.ratio=1,
                           plot.density=TRUE)


}


