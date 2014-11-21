gti<-function(data, annotation=NULL, annID,annName,save.path=NULL,diff=TRUE, cpmtype="Mann-Whitney",
              num.id=TRUE,under=FALSE,p=NULL){

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
  ### GTI Stat ###
  ################
  num = dim(data)[1]
  
  gti = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/num)*(mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))+((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))
  gti_under = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/num)*(mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))-((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))
  
  #Convert data frame into a matrix
  data1<-data.matrix(data)
  
  if(under==FALSE){data2 <- apply_pb(data1,2,gti,Title="Calculating GTI ...")}
  if(under==TRUE){data2 <- apply_pb(data1,2,gti_under,Title="Calculating GTI ...")}
  #data2<-subset(data2,is.nan(data2)==FALSE)
  #data2<-subset(data2,is.na(data2)==FALSE)
  #data2<-subset(data2,data2!="Inf")
  
if(is.null(cpmtype) & !is.null(p)){
    if(under==FALSE){data3<-ifelse(data2>quantile(data2, p,na.rm=TRUE), data2, NA)}
    if(under==TRUE){data3<-ifelse(data2<quantile(data2, p,na.rm=TRUE), data2, NA)}
}
if(!is.null(cpmtype) & is.null(p)){  
  if(under==FALSE){data2.s<-sort(data2,decreasing=TRUE)}
  if(under==TRUE){data2.s<-sort(data2)}
  
  if (diff==FALSE){
    cpm_g<-detectChangePoint(data2.s,cpmType=cpmtype, ARL0=500, startup=20)
  }
  if (diff==TRUE){
    diffg<-c(rep(0,length(data2)))
    for (i in 1:length(data2)-1){
      diffg[i]<-abs(data2.s[i]-data2.s[i+1])
    }
    diffg1<-subset(diffg,is.nan(diffg)==FALSE)
    diffg2<-subset(diffg1,is.na(diffg1)==FALSE)
    cpm_g<-detectChangePoint(diffg2,cpmType="Mann-Whitney", ARL0=500, startup=20)
  }  
  if(under==FALSE){data3<-ifelse(data2>=data2.s[cpm_g$changePoint], data2, NA)}
  if(under==TRUE){data3<-ifelse(data2<=data2.s[cpm_g$changePoint], data2, NA)}
}
data4<-subset(data3, data3!="NA")
  
  if(num.id==TRUE){
    noX<-substr(names(data4), 2,15)
    noX1<-unique(noX)
    if(!is.null(annotation)){
      gsg<-annotation[annotation[,annID] %in% noX1,c(annID,annName)]
      Val<-data4[noX1 %in% gsg[,1]]
      val<-Val[unique(names(Val))]
    }
    else {gsg<-noX1
          val<-as.vector(data4)}
  }
  if(num.id==FALSE){
    #probes<-unique(names(data4))
    if(!is.null(annotation)){
      gsg<-annotation[annotation[,annID] %in% names(data4),c(annID,annName)]
      val<-data4[names(data4) %in% gsg[,1]]
    }
    else {gsg<-names(data4)
          val<-as.vector(data4)}
  }
  gti_out<-data.frame(gsg,Value=val,Stat="GTI")
  rownames(gti_out)<-NULL
  
  if(!is.null(save.path)){
    write.csv(gti_out,file=paste(save.path,"gti_outlier_1grp_cpm.csv"))
  }

## gti cpm plot
  if(!is.null(save.path)){
    setwd(save.path)
    jpeg(filename="gti_plot.jpg")
    
    par(mfrow=c(1,2))
    cpm_plot_D_g<-plot(cpm_g$Ds)
    cpm_thresh_plot_g<-plot(sort(data2),ylab="GTI Values")
    abline(h=data2.s[cpm_g$changePoint],lty=2)
    dev.off()
    #reset window
    par(mfrow=c(1,1))
  }  
  if(is.null(save.path)){
    cpm_thresh_plot_g<-plot(sort(data2),ylab="GTI Values",main="GTI")
    abline(h=data2.s[cpm_g$changePoint],lty=2)
  }
  return(gti_out)
  
}
