
#' Plot the heatmap of the individual PFAS concentrations and PFAS burden scores 

HeatmapPFAS<-function(plot.object){
    
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required but not installed.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required but not installed.")
  }
  if(!is.list(plot.object)){
    stop("Input must be an output object from the PFAScalc2017 function")
  }
  if(is.null(plot.object$use.isomers)){
    stop("Input must be an output object from the PFAScalc2017 function")
  }
  
  plot.data<-plot.object$pfas.burden
  nm.8pfas<-c("me.PFOSA.AcOH","PFUA","PFDeA","PFHxS","PFNA","n.PFOA","Sm.PFOS","n.PFOS")
  nm.7pfas<-c("me.PFOSA.AcOH","PFUA","PFDeA","PFHxS","PFNA","PFOA","PFOS")
  
  if (identical(colnames(plot.data),c(nm.7pfas,"PFAS burden","SE"))){
    ### for 7 PFASs (summed isomers)
    
    plot.data$Summed.concentration<-apply(plot.data[,nm.7pfas],1,function(x){sum(x,na.rm=TRUE)})
    dat.hm<-as.matrix(plot.data[,c(nm.7pfas,"PFAS burden","Summed.concentration")])
    dat.hm<-dat.hm[rev(order(dat.hm[,"PFAS burden"])),]
    
    plot.output<-ComplexHeatmap::Heatmap(dat.hm[,nm.7pfas],
            cluster_columns=F,cluster_rows=F,
            col=circlize::colorRamp2(c(0,10),c("white","darkorange")),
            heatmap_legend_param=list(title="Individual PFAS (ng/mL)"))+
      ComplexHeatmap::Heatmap(dat.hm[,"Summed.concentration"],
              cluster_columns=F,cluster_rows=F,
              col=circlize::colorRamp2(c(0,20),c("white","tan4")),
              heatmap_legend_param=list(title="Summed PFAS (ng/mL)"),
              name="Summed PFAS")+
      ComplexHeatmap::Heatmap(dat.hm[,"PFAS burden"],
              cluster_columns=F,cluster_rows=F,
              col=circlize::colorRamp2(c(-2,0,2),c("blue","white","red")),
              #col=circlize::colorRamp2(c(0,10),c("white","red")),
              heatmap_legend_param=list(title="PFAS exposure burden"),
              name="PFAS burden")
    
  }
  else{
    ### for 8 PFASs (isomers)
    plot.data$Summed.concentration<-apply(plot.data[,nm.8pfas],1,function(x){sum(x,na.rm=TRUE)})
    dat.hm<-as.matrix(plot.data[,c(nm.8pfas,"PFAS burden (isomers)","Summed.concentration")])
    dat.hm<-dat.hm[rev(order(dat.hm[,"PFAS burden (isomers)"])),]
    
    plot.output<-ComplexHeatmap::Heatmap(dat.hm[,nm.8pfas],
                         cluster_columns=F,cluster_rows=F,
                         col=circlize::colorRamp2(c(0,10),c("white","darkorange")),
                         heatmap_legend_param=list(title="Individual PFAS (ng/mL)"))+
      ComplexHeatmap::Heatmap(dat.hm[,"Summed.concentration"],
              cluster_columns=F,cluster_rows=F,
              col=circlize::colorRamp2(c(0,20),c("white","tan4")),
              heatmap_legend_param=list(title="Summed PFAS (ng/mL)"),
              name="Summed PFAS")+
      ComplexHeatmap::Heatmap(dat.hm[,"PFAS burden (isomers)"],
              cluster_columns=F,cluster_rows=F,
              col=circlize::colorRamp2(c(-2,0,2),c("blue","white","red")),
              #col=circlize::colorRamp2(c(0,10),c("white","red")),
              heatmap_legend_param=list(title="PFAS exposure burden \n (isomers)"),
              name="PFAS burden \n (isomers)")
    
  }
  
  return(plot.output)
  
}
