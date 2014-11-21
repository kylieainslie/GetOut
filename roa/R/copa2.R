copa2<-function(data_d,data_n,data_all, annotation=NULL, annID,annName,save.path=NULL,diff=TRUE, cpmtype="Mann-Whitney",
                cut=0.85,num.id=TRUE,under=FALSE,p=NULL){

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
  
  
  #check if cut-off matches under status
  if(under==TRUE & cut>0.5){warning("cut value doesn't match under status")}
  if(under==FALSE & cut<0.5){warning("cut value doesn't match under status")}
  
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
#################
### COPA Stat ###
#################

#COPA             
# use combined (normal and disease) samples for median and mad

#take rth percentile of disease group
copa<-function(data_d,data_n,data_all,Cut=cut){
  #calculate quantile specified by cut for each gene
  r<-function(data) (quantile(data,Cut,na.rm=TRUE)) #quantile function
  x<-data.matrix(data_d)
  x1<-apply_pb(x,2,r)
  #calculate median
  #median1<- function(data) quantile(data,0.5,na.rm=TRUE)
  data1<-data.matrix(data_all)
  med<-apply_pb(data1,2,median,Title="Calculating medians ...")
  #calculate mad
  #mad1<-function(data) mad(data,na.rm=TRUE)
  mad1<-apply_pb(data1,2,mad,Title="Calculating median absolute deviations ...")
  #calculate copa
  copa<-(x1-med)/mad1
  return(copa)
}

  #Convert data frame into a matrix
  data1<-data.matrix(data_d)
  data2<-data.matrix(data_n)
  data3<-data.matrix(data_all)

  
  data2c <- copa(data1,data2,data3,Cut=cut)
  data2c<-subset(data2c,is.nan(data2c)==FALSE)
  data2c<-subset(data2c,is.na(data2c)==FALSE)
  data2c<-subset(data2c,data2c!="Inf")

if(is.null(cpmtype) & !is.null(p)){
  if(under==FALSE){data3<-ifelse(data2c>quantile(data2c, p,na.rm=TRUE), data2c, NA)}
  if(under==TRUE){data3<-ifelse(data2c<quantile(data2c, p,na.rm=TRUE), data2c, NA)}
}
if(!is.null(cpmtype) & is.null(p)){
  if(under==FALSE){data2c.s<-sort(data2c,decreasing=TRUE)}
  if(under==TRUE){data2c.s<-sort(data2c)}
  
  if (diff==FALSE){
    cpm_c<-detectChangePoint(data2c.s,cpmType=cpmtype, ARL0=500, startup=20)
  }
  if (diff==TRUE){
    diffc<-c(rep(0,length(data2c)))
    for (i in 1:length(data2c)-1){
      diffc[i]<-abs(data2c.s[i]-data2c.s[i+1])
    }
    diffc<-subset(data2c,is.nan(diffc)==FALSE)
    diffc<-subset(data2c,is.na(diffc)==FALSE)
    cpm_c<-detectChangePoint(diffc,cpmType="Mann-Whitney", ARL0=500, startup=20)
  }  
  if(under==FALSE){data3c<-ifelse(data2c>=data2c.s[cpm_c$changePoint], data2c, NA)}
  if(under==TRUE){data3c<-ifelse(data2c<=data2c.s[cpm_c$changePoint], data2c, NA)}
}
data4c<-subset(data3c, data3c!="NA")
  
if(num.id==TRUE){
  noX<-substr(names(data4c), 2,15)
  noX1<-unique(noX)
  if(!is.null(annotation)){
    gsg<-annotation[annotation[,annID] %in% noX1,c(annID,annName)]
    Val<-data4c[noX1 %in% gsg[,1]]
    val<-Val[unique(names(Val))]
  }
  else {gsg<-noX1
        val<-as.vector(data4c)}
}
if(num.id==FALSE){
  if(!is.null(annotation)){
    gsg<-annotation[annotation[,annID] %in% names(data4c),c(annID,annName)]
    val<-data4c[names(data4c) %in% gsg[,1]]
  }
  else {gsg<-names(data4c)
        val<-as.vector(data4c)}
}
copa_out<-data.frame(gsg,Value=val,Stat="COPA")
rownames(copa_out)<-NULL

if(!is.null(save.path)){
  write.csv(copa_out,file=paste(save.path,"copa_outlier_1grp_cpm.csv"))
}

## COPA cpm plot
if(!is.null(save.path)){
  setwd(save.path)
  jpeg(filename="copa_plot.jpg")
  
  par(mfrow=c(1,2))
  cpm_plot_D_c<-plot(cpm_c$Ds)
  cpm_thresh_plot_c<-plot(sort(data2c),ylab="COPA Values")
  abline(h=data2c.s[cpm_c$changePoint],lty=2)
  dev.off()
  #reset window
  par(mfrow=c(1,1))
}
if(is.null(save.path)){
  cpm_thresh_plot_c<-plot(sort(data2c),ylab="COPA Values",main="COPA")
  abline(h=data2c.s[cpm_c$changePoint],lty=2) 
}
## Output
  return(copa_out)
  
}


