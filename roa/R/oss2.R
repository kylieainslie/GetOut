oss2<-function(data_d,data_n,data_all, annotation=NULL, annID,annName,save.path=NULL,diff=T, cpmtype="Mann-Whitney",
               num.id=T,under=F,p=NULL){
  
  #check for missing arguments
  if(missing(data_d) | missing(data_n) | missing(data_all)){
    stop("No data specified.")
  }
  #check for numeric arguments
  if(!is.numeric(data_d) | !is.numeric(data_n) | !is.numeric(data_all)){
    stop("Need numeric data matrix.")
  }
  #check that columns for annID and annName are numeric and not larger than total number of columns
  if(is.null(annotation)){
    annID=NULL
    annName=NULL
  }
  else if((!is.null(annotation) & !is.numeric(annID)) | (!is.null(annotation) & !is.numeric(annName))
          | (!is.null(annotation) & annID>dim(annotation)[2]) | (!is.null(annotation) & annName>dim(annotation)[2])){
    stop("Invalid annotation column specification")
  }  
  
#both cpmtype and p cannot be null
  if(!is.null(cpmtype)){p=NULL}
  if(!is.null(p)){cpmtype=NULL}
  if(is.null(cpmtype) & is.null(p)){stop("Must select change point model or user defined cut-off.")}
  
  ####################
  ### Progress Bar ###  
  ####################
  apply_pb <- function(X, MARGIN, FUN, ...,Title="Calculating ...")
  {
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- tkProgressBar(title=Title,min=0, max=pb_Total,width=300)
    
    wrapper <- function(...)
    {
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir= env)
      setTkProgressBar(get("pb", envir= env),curVal +1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
###########
### OSS ###
###########
#OSS
oss<-function(data_d,data_n,data_all,Under=under){
  #calculate median
  data1<-data.matrix(data_all)
  med<-apply_pb(data1,2,median,Title="Calculating medians ...")
  #calculate mad
  data2<-data.matrix(data_n)
  mad2<-apply_pb(data2,2,mad,Title="Calculating median absolute deviations ...")
  #standardize gene expression
  x2<-sweep(data_all,2,med,FUN="-")
  x3<-sweep(x2,2,mad2,FUN="/")
  #determine cutoff
  if(Under==F){cutoff<-function(data) quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T)+quantile(data,0.75,na.rm=T)}
  if(Under==T){cutoff<-function(data) quantile(data,0.25,na.rm=T)-(quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))}
  
  gene_cut<-apply_pb(data.matrix(x3),2,cutoff,Title="Calculating cut-offs ...")
  
  #standardize disease group gene expression
  y2<-sweep(data_d,2,med,FUN="-")
  y3<-sweep(x2,2,mad2,FUN="/")
  y4<-data.matrix(y3)   
  #calculate oss
  oss<-c(rep(0,dim(y4)[2]))
  for (i in 1:dim(y4)[2]){
    if(Under==F){oss[i]<-sum(subset(y4[,i],y4[,i]>gene_cut[i]))}
    if(Under==T){oss[i]<-sum(subset(y4[,i],y4[,i]<gene_cut[i]))}
  }
  names(oss)<-colnames(data_all)
  return(oss)
}

################################
### Apply Outlier Statistics ###
################################
data1<-data.matrix(data_d)
data2<-data.matrix(data_n)
data3<-data.matrix(data_all)

data2o<-oss(data1,data2,data3)

#####################################
### Change point model and output ###
#####################################
if(is.null(cpmtype) & !is.null(p)){
  if(under==F){data3o<-ifelse(data2o>quantile(data2o, p,na.rm=TRUE), data2o, NA)}
  if(under==T){data3o<-ifelse(data2o<quantile(data2o, p,na.rm=TRUE), data2o, NA)}
}
if(!is.null(cpmtype) & is.null(p)){  
  if(under==F){data2o.s<-sort(data2o,decreasing=T)}
  if(under==T){data2o.s<-sort(data2o)}
  
  if (diff==F){
    cpm_o<-detectChangePoint(data2o.s,cpmType="Mann-Whitney", ARL0=500, startup=6)
  }
  if (diff==T){
    diffo<-c(rep(0,length(data2o)))
    for (i in 1:length(data2o)-1){
      diffo[i]<-abs(data2o.s[i]-data2o.s[i+1])
    }
    diffo<-subset(data2o,is.nan(diffo)==F)
    diffo<-subset(data2o,is.na(diffo)==F)
    cpm_o<-detectChangePoint(diffo,cpmType="Mann-Whitney", ARL0=500, startup=6)
  }  
  if(under==F){data3o<-ifelse(data2o>=data2o.s[cpm_o$changePoint], data2o, NA)}
  if(under==T){data3o<-ifelse(data2o<=data2o.s[cpm_o$changePoint], data2o, NA)}
}  
data4o<-subset(data3o, data3o!="NA")

if(num.id==T){
  noX<-substr(names(data4o), 2,15)
  noX1<-unique(noX)
  if(!is.null(annotation)){
    gso<-annotation[annotation[,annID] %in% noX1,c(annID,annName)]
    Val<-data4o[noX1 %in% gso[,1]]
    val<-Val[unique(names(Val))]
  }
  else {gso<-noX1
        val<-as.vector(data4o)}
}
if(num.id==F){
  #probes<-unique(names(data4))
  if(!is.null(annotation)){
    gso<-annotation[annotation[,annID] %in% names(data4o),c(annID,annName)]
    val<-data4o[names(data4o) %in% gso[,1]]
  }
  else {gso<-names(data4o)
        val<-as.vector(data4o)}
}
oss_out<-data.frame(gso,Value=val,Stat="OSS")
rownames(oss_out)<-NULL

if(!is.null(save.path)){
  write.csv(oss_out,file=paste(save.path,"oss_outlier_1grp_cpm.csv"))
}

##OSS cpm plot
if(!is.null(save.path)){
  setwd(save.path)
  jpeg(filename="oss_plot.jpg")
  
  par(mfrow=c(1,2))
  cpm_plot_D_o<-plot(cpm_o$Ds)
  cpm_thresh_plot_o<-plot(sort(data2o),ylab="OSS Values")
  abline(h=data2o.s[cpm_o$changePoint],lty=2)
  dev.off()
  #reset window
  par(mfrow=c(1,1))
} 
if(is.null(save.path)){
  cpm_thresh_plot_o<-plot(sort(data2o),ylab="OSS Values",main="OSS")
  abline(h=data2o.s[cpm_o$changePoint],lty=2)
}
return(oss_out)

}
