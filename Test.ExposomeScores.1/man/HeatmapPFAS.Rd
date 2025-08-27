\name{HeatmapPFAS}
\alias{HeatmapPFAS}
\title{Plotting the heatmap of the individual PFAS concentrations and PFAS burden scores}
\usage{
HeatmapPFAS(plot.object)
}
\description{
Plotting the heatmap of the individual PFAS concentrations and PFAS burden scores.
}

\arguments{
  \item{plot.object}{An output object from the PFAScalc2017 function.}
}


\value{
A Heatmap-class object.
}
\examples{
PFAS.burden<-PFAScalc2017(data=PFAScycleJ,
             me.PFOSA.AcOH="me.PFOSA.AcOH",
             PFUA="PFUA",
             PFDeA="PFDeA",
             PFHxS="PFHxS",
             PFNA="PFNA",
             PFOA="PFOA",
             PFOS="PFOS",
             use.isomers=FALSE)

HeatmapPFAS(PFAS.burden)

}


