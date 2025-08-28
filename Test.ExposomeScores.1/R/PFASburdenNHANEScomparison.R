
#' Comparing the estimated PFAS burden score with 2017-2018 NHANES

PFASburdenNHANEScomparison<-function(ID,Sex,Age,
                                          PFAS.burden,use.isomers=FALSE,
                                          matched.ratio=5,
                                          bins.num=50,
                                          plot.density=TRUE){
    
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  
  #require(MatchIt)
  require(ggplot2)
  
  data("PFAScycleJ", envir = environment())
  
  data.1<-data.frame(ID=ID,Male=Sex,Age=Age,PFAS.burden=PFAS.burden)
  data.1$Sample<-1
  
  if (use.isomers){
    ### for 8 PFASs (using isomers)
    data.0<-data.frame(ID=PFAScycleJ$SEQN,
                       Male=PFAScycleJ$Male,
                       Age=PFAScycleJ$Age,
                       PFAS.burden=PFAScycleJ$`PFAS burden (isomers)`)
    data.0$Sample<-0
    
    data.01<-rbind(data.0,data.1)
    mtch.01<-MatchIt::matchit(Sample~Male+Age,data=data.01,
                              exact=~Male,
                              distance="mahalanobis",method="nearest",
                              ratio=matched.ratio,replace=TRUE)
    mtch.data<-MatchIt::get_matches(mtch.01)
    mtch.data$Comparison<-factor(dplyr::recode(mtch.data$Sample,
                                               "0"="2017-2018 NHANES",
                                               "1"="Study Sample"),
                                 levels=c("2017-2018 NHANES","Study Sample"))
    
    output.histogram<-ggplot(mtch.data,aes(x=PFAS.burden,color=Comparison,fill=Comparison))+
      geom_histogram(alpha=0.5,position="identity",bins=bins.num)+
      scale_fill_manual(values=c("dodgerblue","coral"))+
      scale_color_manual(values=c("dodgerblue","coral"))+
      labs(x="Estimated PFAS burden (isomers)",
           title="Comparing with the 2017-2018 NHANES \nafter matching for sex and age")+
      theme_bw()
    
    output.density<-ggplot(mtch.data,aes(x=PFAS.burden,color=Comparison,fill=Comparison))+
      geom_density(alpha=0.3)+
      scale_fill_manual(values=c("dodgerblue","coral"))+
      scale_color_manual(values=c("dodgerblue","coral"))+
      labs(x="Estimated PFAS burden (isomers)",
           title="Comparing with the 2017-2018 NHANES \nafter matching for sex and age")+
      theme_bw()
    
  }
  else{
    ### for 7 PFASs (using sum of isomers)
    data.0<-data.frame(ID=PFAScycleJ$SEQN,
                       Male=PFAScycleJ$Male,
                       Age=PFAScycleJ$Age,
                       PFAS.burden=PFAScycleJ$`PFAS burden`)
    data.0$Sample<-0
    
    data.01<-rbind(data.0,data.1)
    mtch.01<-MatchIt::matchit(Sample~Male+Age,data=data.01,
                              exact=~Male,
                              distance="mahalanobis",method="nearest",
                              ratio=matched.ratio,replace=TRUE)
    mtch.data<-MatchIt::get_matches(mtch.01)
    mtch.data$Comparison<-factor(dplyr::recode(mtch.data$Sample,
                                               "0"="2017-2018 NHANES",
                                               "1"="Study Sample"),
                                 levels=c("2017-2018 NHANES","Study Sample"))
    
    output.histogram<-ggplot(mtch.data,aes(x=PFAS.burden,color=Comparison,fill=Comparison))+
      geom_histogram(alpha=0.5,position="identity",bins=bins.num)+
      scale_fill_manual(values=c("dodgerblue","coral"))+
      scale_color_manual(values=c("dodgerblue","coral"))+
      labs(x="Estimated PFAS burden",
           title="Comparing with the 2017-2018 NHANES \nafter matching for sex and age")+
      theme_bw()
    
    output.density<-ggplot(mtch.data,aes(x=PFAS.burden,color=Comparison,fill=Comparison))+
      geom_density(alpha=0.3)+
      scale_fill_manual(values=c("dodgerblue","coral"))+
      scale_color_manual(values=c("dodgerblue","coral"))+
      labs(x="Estimated PFAS burden",
           title="Comparing with the 2017-2018 NHANES \nafter matching for sex and age")+
      theme_bw()
    
    
    
  }
  
  if(plot.density){
    plot.output<-output.density
  } else {
    plot.output<-output.histogram
  }
  
  return(plot.output)
  
}
