PFAScalc2017<-function(data,me.PFOSA.AcOH,PFUA,PFDeA,PFHxS,PFNA,PFOA,PFOS){
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' is required but not installed.")
  }
  
  data("WtdQuantile", package = "Test.ExposomeScores", envir = environment())
  data("mod.mirt.sum", package = "Test.ExposomeScores", envir = environment())
  data("mod.mirt.isomer", package = "Test.ExposomeScores", envir = environment())
  
  data.cont<-data[,c(me.PFOSA.AcOH,PFUA,PFDeA,PFHxS,PFNA,PFOA,PFOS)]
  nm.7pfas<-c("me.PFOSA.AcOH","PFUA","PFDeA","PFHxS","PFNA","PFOA","PFOS")
  colnames(data.cont)<-nm.7pfas
  
  data.cut<-data.cont
  
  for (pfas.i in nm.7pfas){
    data.cut[,pfas.i]=as.numeric(cut(data.cont[,pfas.i],WtdQuantile[[pfas.i]],right=FALSE,labels=c(1:(length(WtdQuantile[[pfas.i]])-1))))
  }
  
  
  fs.output<-mirt::fscores(mod.mirt.sum,response.pattern=data.cut)
  return(fs.output)
  
  
}