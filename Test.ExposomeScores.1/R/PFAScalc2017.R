#' Calculate PFAS exposure burden
#'
#' This function calculate PFAS exposure burden using the 2017 - 2018 US PFAS exposure burden calculator.
#' 
#' @param data A dataframe of individual PFAS measurements
#' @param me.PFOSA.AcOH The column name of me-PFOSA-AcOH in the dataframe.
#' @param PFUA The column name of PFUA in the dataframe.
#' @param PFDeA The column name of PFDeA in the dataframe.
#' @param PFHxS The column name of PFHxS in the dataframe.
#' @param PFNA The column name of PFNA in the dataframe.
#' @param PFOA The column name of PFOA in the dataframe.
#' @param PFOS The column name of PFOS in the dataframe.
#' @param n.PFOA The column name of n-PFOA in the dataframe.
#' @param Sm.PFOS The column name of Sm-PFOS in the dataframe.
#' @param n.PFOS The column name of n-PFOS in the dataframe.
#' @param use.isomers Whether using isomers for PFOA and PFOS
#' 
#' @return A dataframe with individual PFAS measurements and estimated PFAS exposure burden and SE of the estimated PFAS burden

PFAScalc2017<-function(data,
                       me.PFOSA.AcOH=NA,PFUA=NA,PFDeA=NA,PFHxS=NA,PFNA=NA,
                       PFOA=NA,PFOS=NA,
                       n.PFOA=NA,Sm.PFOS=NA,n.PFOS=NA,
                       use.isomers=FALSE){
  
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' is required but not installed.")
  }
  
  data("WtdQuantile", envir = environment())
  data("mod.mirt.sum",envir = environment())
  data("mod.mirt.isomer",envir = environment())
  
    
  if (use.isomers) {
    nm.8pfas<-c("me.PFOSA.AcOH","PFUA","PFDeA","PFHxS","PFNA","n.PFOA","Sm.PFOS","n.PFOS")
    
    if (!is.na(me.PFOSA.AcOH)){conc.me.PFOSA.AcOH<-data[,me.PFOSA.AcOH]} else {conc.me.PFOSA.AcOH<-rep(NA,nrow(data))}
    if (!is.na(PFUA)){conc.PFUA<-data[,PFUA]} else {conc.PFUA<-rep(NA,nrow(data))}
    if (!is.na(PFDeA)){conc.PFDeA<-data[,PFDeA]} else {conc.PFDeA<-rep(NA,nrow(data))}
    if (!is.na(PFHxS)){conc.PFHxS<-data[,PFHxS]} else {conc.PFHxS<-rep(NA,nrow(data))}
    if (!is.na(PFNA)){conc.PFNA<-data[,PFNA]} else {conc.PFNA<-rep(NA,nrow(data))}
    if (!is.na(n.PFOA)){conc.n.PFOA<-data[,n.PFOA]} else {conc.n.PFOA<-rep(NA,nrow(data))}
    if (!is.na(Sm.PFOS)){conc.Sm.PFOS<-data[,Sm.PFOS]} else {conc.Sm.PFOS<-rep(NA,nrow(data))}
    if (!is.na(n.PFOS)){conc.n.PFOS<-data[,n.PFOS]} else {conc.n.PFOS<-rep(NA,nrow(data))}
    
    data.cont<-data.frame(me.PFOSA.AcOH=conc.me.PFOSA.AcOH,
                          PFUA=conc.PFUA,
                          PFDeA=conc.PFDeA,
                          PFHxS=conc.PFHxS,
                          PFNA=conc.PFNA,
                          n.PFOA=conc.n.PFOA,
                          Sm.PFOS=conc.Sm.PFOS,
                          n.PFOS=conc.n.PFOS)
    
    data.cut<-data.cont
    for (pfas.i in nm.8pfas){
      try((data.cut[,pfas.i]=as.numeric(cut(data.cont[,pfas.i],WtdQuantile[[pfas.i]],right=FALSE,labels=c(1:(length(WtdQuantile[[pfas.i]])-1))))),silent=TRUE)
    }
    
    fs<-mirt::fscores(mod.mirt.isomer,response.pattern=data.cut)
    colnames(fs)<-c("PFAS burden (isomers)","SE")
    fs.output<-cbind(data.cont,fs)

  }
  else {
    nm.7pfas<-c("me.PFOSA.AcOH","PFUA","PFDeA","PFHxS","PFNA","PFOA","PFOS")
    
    if (!is.na(me.PFOSA.AcOH)){conc.me.PFOSA.AcOH<-data[,me.PFOSA.AcOH]} else {conc.me.PFOSA.AcOH<-rep(NA,nrow(data))}
    if (!is.na(PFUA)){conc.PFUA<-data[,PFUA]} else {conc.PFUA<-rep(NA,nrow(data))}
    if (!is.na(PFDeA)){conc.PFDeA<-data[,PFDeA]} else {conc.PFDeA<-rep(NA,nrow(data))}
    if (!is.na(PFHxS)){conc.PFHxS<-data[,PFHxS]} else {conc.PFHxS<-rep(NA,nrow(data))}
    if (!is.na(PFNA)){conc.PFNA<-data[,PFNA]} else {conc.PFNA<-rep(NA,nrow(data))}
    if (!is.na(PFOA)){conc.PFOA<-data[,PFOA]} else {conc.PFOA<-rep(NA,nrow(data))}
    if (!is.na(PFOS)){conc.PFOS<-data[,PFOS]} else {conc.PFOS<-rep(NA,nrow(data))}
    
    data.cont<-data.frame(me.PFOSA.AcOH=conc.me.PFOSA.AcOH,
                          PFUA=conc.PFUA,
                          PFDeA=conc.PFDeA,
                          PFHxS=conc.PFHxS,
                          PFNA=conc.PFNA,
                          PFOA=conc.PFOA,
                          PFOS=conc.PFOS)
    data.cut<-data.cont
    for (pfas.i in nm.7pfas){
      try((data.cut[,pfas.i]=as.numeric(cut(data.cont[,pfas.i],WtdQuantile[[pfas.i]],right=FALSE,labels=c(1:(length(WtdQuantile[[pfas.i]])-1))))),silent=TRUE)
    }
    
    fs<-mirt::fscores(mod.mirt.sum,response.pattern=data.cut)
    colnames(fs)<-c("PFAS burden","SE")
    fs.output<-cbind(data.cont,fs)
    
  }
  return(fs.output)
}