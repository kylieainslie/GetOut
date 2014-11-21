variance<-function(data, annotation=NULL, annID,annName,save.path=NULL,diff=TRUE, cpmtype="Mann-Whitney",cut=0.85,
                  num.id=TRUE,p=NULL){

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
  ### Var Stat ###
  ################

  var.stat = function(data) var(data)

  #Convert data frame into a matrix
  data1<-data.matrix(data)

  data2v <- apply_pb(data1,2,var.stat,Title="Calculating variance ...")
  data2v<-subset(data2v,is.nan(data2v)==FALSE)
  data2v.s<-sort(data2v,decreasing=TRUE)
  
if(is.null(cpmtype) & !is.null(p)){
    data3v<-ifelse(data2v>quantile(data2v, p,na.rm=TRUE), data2v, NA)
}
  
if(!is.null(cpmtype) & is.null(p)){
  if (diff==FALSE){
    cpm_v<-detectChangePoint(data2v.s,cpmType="Mann-Whitney", ARL0=500, startup=20)
  }
  if (diff==TRUE){
  diffv<-c(rep(0,length(data2v)))
    for (i in 1:length(data2v)-1){
    diffv[i]<-abs(data2v.s[i]-data2v.s[i+1])
    }
  cpm_v<-detectChangePoint(diffv,cpmType="Mann-Whitney", ARL0=500, startup=20)
  }  
  data3v<-ifelse(data2v>=data2v.s[cpm_v$changePoint], data2v, NA)
}
data4v<-subset(data3v, data3v!="NA")

  if(num.id==TRUE){
    noX<-substr(names(data4v), 2,15)
    noX1<-unique(noX)
    if(!is.null(annotation)){
      gsv<-annotation[annotation[,annID] %in% noX1,c(annID,annName)]
      Val<-data4v[noX1 %in% gsv[,1]]
      val<-Val[unique(names(Val))]
    }
    else {gsv<-noX1
          val<-as.vector(data4v)}
  }
  if(num.id==FALSE){
    #probes<-unique(names(data4))
    if(!is.null(annotation)){
      gsv<-annotation[annotation[,annID] %in% names(data4v),c(annID,annName)]
      val<-data4v[names(data4v) %in% gsv[,1]]
    }
    else {gsv<-names(data4v)
          val<-as.vector(data4v)}
  }
  var_out<-data.frame(gsv,Value=val,Stat="Var")
  rownames(var_out)<-NULL
  
  if(!is.null(save.path)){
    write.csv(var_out,file=paste(save.path,"var_outlier_1grp_cpm.csv"))
  }
#  if (outsample==TRUE & !is.null(save.path)){
#    outlying.samples(data,names(data4),over_cut=cut,stat="Var",savepath=save.path) 
#  }

##Var cpm plot
if(!is.null(save.path)){
  setwd(save.path)
  jpeg(filename="var_plot.jpg")
  
  par(mfrow=c(1,2))
  cpm_plot_D_v<-plot(cpm_v$Ds)
  cpm_thresh_plot_v<-plot(sort(data2v),ylab="Variance Values",main="Variance")
  abline(h=data2v.s[cpm_v$changePoint],lty=2)
  dev.off()
  #reset window
  par(mfrow=c(1,1))
}  
if(is.null(save.path)){
  cpm_thresh_plot_v<-plot(sort(data2v),ylab="Variance Values",main="Variance")
  abline(h=data2v.s[cpm_v$changePoint],lty=2) 
}
return(var_out)

}
