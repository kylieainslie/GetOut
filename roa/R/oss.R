oss<-function(data, annotation=NULL, annID,annName,save.path=NULL,diff=T, cpmtype="Mann-Whitney",
              num.id=T,under=F,p=NULL){

  #check for missing arguments
  if(missing(data)){
    stop("No data specified.")
  }
  #check for numeric arguments
  if(!is.numeric(data)){
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
  ################
  ### OSS Stat ###
  ################
  
  oss = function(data) sum(subset(data,data>((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T)))-median(data,na.rm=T))/mad(data,na.rm=T)
  oss_under = function(data) sum(subset(data,data<(quantile(data,0.25,na.rm=T)-(quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))))-median(data))/mad(data)
  
  #Convert data frame into a matrix
  data1<-data.matrix(data)
  
  if(under==F){data2o <- apply_pb(data1,2,oss,Title="Calculating OSS ...")}
  if(under==T){data2o <- apply_pb(data1,2,oss_under,Title="Calculating OSS ...")}
  #data2o<-subset(data2o,is.nan(data2o)==F)
  #data2o<-subset(data2o,is.na(data2o)==F)
  #data2o<-subset(data2o,data2o!="Inf")
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
