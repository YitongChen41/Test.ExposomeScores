
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
    
    fs.output<-mirt::fscores(mod.mirt.isomer,response.pattern=data.cut)
    colnames(fs.output)<-c("PFAS burden (isomers)","SE")
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
    
    fs.output<-mirt::fscores(mod.mirt.sum,response.pattern=data.cut)
    colnames(fs.output)<-c("PFAS burden","SE")
    
  }
  return(fs.output)
}