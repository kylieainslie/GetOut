gti2<-function(data_d,data_n,data_all, annotation=NULL, annID,annName,save.path=NULL,diff=TRUE, 
               cpmtype="Mann-Whitney",num.id=TRUE,under=FALSE,p=NULL){

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
  ################
  ### GTI Stat ###
  ################
  ##Define number of samples per probe
  
  num1 = dim(data_n)[1]
  num2 = dim(data_d)[1]
  
  #GTI - a statistic is calculated for each group (1=normal; 2=disease)
 gti1 = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/num1)*(mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))+((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))  
 gti2 = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/num2)*(mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))+((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))  
 
 gti_under1 = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/num1)*(mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))-((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))
 gti_under2 = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/num2)*(mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))-((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))
 
  gti<-function(data_d,data_n,data_all,Under=under){
    if(Under==FALSE){
      data_norm <- apply_pb(data1,2,gti1, Title="Calculating GTI ...")
      data_dis<- apply_pb(data2,2,gti2, Title="Calculating GTI ...")
    }
    if(Under==TRUE){
      data_norm <- apply_pb(data1,2,gti_under1, Title="Calculating GTI ...")
      data_dis<- apply_pb(data2,2,gti_under2, Title="Calculating GTI ...")
    }
    dist<-data_dis-data_norm
    return(dist)  
  }
  
  #Convert data frame into a matrix
  data1<-data.matrix(data_d)
  data2<-data.matrix(data_n)
  data3<-data.matrix(data_all)
  
  data2g<-gti(data1,data2,data3)
  
if(is.null(cpmtype) & !is.null(p)){
    if(under==FALSE){data3<-ifelse(data2g>quantile(data2g, p,na.rm=TRUE), data2g, NA)}
    if(under==TRUE){data3<-ifelse(data2g<quantile(data2g, p,na.rm=TRUE), data2g, NA)}
}
if(!is.null(cpmtype) & is.null(p)){  
  if(under==FALSE){data2g.s<-sort(data2g,decreasing=TRUE)}
  if(under==TRUE){data2g.s<-sort(data2g)}
  
  if (diff==FALSE){
    cpm_g<-detectChangePoint(data2g.s,cpmType=cpmtype, ARL0=500, startup=20)
  }
  if (diff==TRUE){
    diffg<-c(rep(0,length(data2g)))
    for (i in 1:length(data2g)-1){
      diffg[i]<-abs(data2g.s[i]-data2g.s[i+1])
    }
    diffg1<-subset(diffg,is.nan(diffg)==FALSE)
    diffg2<-subset(diffg1,is.na(diffg1)==FALSE)
    cpm_g<-detectChangePoint(diffg2,cpmType="Mann-Whitney", ARL0=500, startup=20)
  }  
  if(under==FALSE){data3<-ifelse(data2g>=data2g.s[cpm_g$changePoint], data2g, NA)}
  if(under==TRUE){data3<-ifelse(data2g<=data2g.s[cpm_g$changePoint], data2g, NA)}
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
  cpm_thresh_plot_g<-plot(sort(data2g),ylab="GTI Values")
  abline(h=data2g.s[cpm_g$changePoint],lty=2)
  dev.off()
  #reset window
  par(mfrow=c(1,1))
}  
if(is.null(save.path)){
  cpm_thresh_plot_g<-plot(sort(data2g),ylab="GTI Values",main="GTI")
  abline(h=data2g.s[cpm_g$changePoint],lty=2)
}
return(gti_out)
  
}
