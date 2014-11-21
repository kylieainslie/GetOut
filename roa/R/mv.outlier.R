######################################
### Multivariate Outlier statistic ###
######################################
mv.outlier<-function(data,annotation=NULL,annID,annName,save.path=NULL,cut=0.85,
                     cpmtype="Mann-Whitney",num.id=TRUE, under=FALSE, diff=TRUE, p=NULL){

  #check for missing arguments
  if(missing(data)){
    stop("No data specified")
  }
  #check for numeric arguments
  if(!is.numeric(data)){
    stop("Need numeric data matrix")
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
  #if(under==TRUE & cut>0.5){warning("cut value doesn't match under status")}
  #if(under==FALSE & cut<0.5){warning("cut value doesn't match under status")}
  
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

#####################
### Outlier Stats ###
#####################
num = dim(data)[1]
##GTI
gti = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/num)*(mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))+((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))
gti_under = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/num)*(mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))-((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))
  
##COPA
copa = function(data) (quantile(data,cut,na.rm=TRUE)-median(data))/mad(data)
##OSS
oss_new = function(data) sum(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))-median(data))/mad(data)
oss_under = function(data) sum(subset(data,data<(quantile(data,0.25,na.rm=TRUE)-(quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))))-median(data))/mad(data)
  
#Convert data frame into a matrix
data1<-data.matrix(data)

#apply outlier stats
if(under==FALSE){data2 <- apply_pb(data1,2,gti, Title="Calculating GTI ...")}
if(under==TRUE){data2 <- apply_pb(data1,2,gti_under, Title="Calculating GTI ...")}
  
data2c <- apply_pb(data1,2,copa, Title="Calculating COPA ...")
  
if(under==FALSE){data2o <- apply_pb(data1,2,oss_new, Title="Calculating OSS ...")}
if(under==TRUE){data2o <- apply_pb(data1,2,oss_under, Title="Calculating OSS ...")}

#calculate score statistic
score<-c(rep(0,length(data2)))
for(i in 1:length(data2)){
  #print(i)
  value<-matrix(c(data2[i],data2c[i],data2o[i]),nrow=3)
  score[i]<-t(value)%*%diag(3)%*%value
}  
score1<-data.matrix(score)
names(score1)<-names(data2)


############################
### Identifying Outliers ###
############################
#determining outliers from cpm
if(!is.null(cpmtype) & is.null(p)){
  score1.s<-sort(score1,decreasing=TRUE)
  
  if (diff==FALSE){cpm_s<-detectChangePoint(score1.s,cpmType=cpmtype, ARL0=500, startup=20)}
  if (diff==TRUE){
    diffs<-c(rep(0,length(score1.s)))
    for (i in 1:length(score1.s)-1){
      diffs[i]<-score1.s[i]-score1.s[i+1]
    }
    diffs1<-subset(diffs,is.nan(diffs)==FALSE)
    diffs2<-subset(diffs1,is.na(diffs1)==FALSE)
    diffs3<-subset(diffs2,diffs2!="Inf")
    
  cpm_s<-detectChangePoint(diffs3,cpmType=cpmtype, ARL0=500, startup=20)
  }
  position<-which(score1>=score1.s[cpm_s$changePoint])
  outliers<-score1[position]
}
  
#determining outliers from user-defined cut-off  
if(is.null(cpmtype) & !is.null(p)){
  score2<-ifelse(score1>quantile(score, p,na.rm=TRUE), score1, NA)
  outliers<-score2[score2!="NA"]
}
##############
### Output ###
##############
  
if(num.id==TRUE){
id<-substr(names(outliers),2,15)
if(is.null(annotation)){
  gene.symbol<-id
}
  else{gene.symbol<-annotation[annotation[,annID]%in% substr(names(outliers),2,15),annName]}
}
if(num.id==FALSE){
id<-names(outliers)
if(is.null(annotation)){
  gene.symbol<-id
}
  else{gene.symbol<-annotation[annotation[,annID]%in% names(outliers),annName]}
}
out.data<-data.frame(ID=id,Gene=gene.symbol,Score=outliers)
rownames(out.data)<-NULL
if(!is.null(save.path)){
setwd(save.path)
write.csv(out.data,paste("mv_outliers_",cut,".csv",sep=""))
}
#############
### Plots ###
#############
if(!is.null(save.path)){
  setwd(save.path)
  jpeg(filename="mv_outlier_plot.jpg")
  
  #par(mfrow=c(1,2))
  cpm_plot_D_v<-plot(cpm_s$Ds)
  cpm_thresh_plot_s<-plot(score1.s,ylab="Score Values",main="Multivariate Score")
  abline(h=score1.s[cpm_s$changePoint],lty=2)
  dev.off()
  #reset window
  #par(mfrow=c(1,1))
}  
if(is.null(save.path)){
  cpm_thresh_plot_s<-plot(score1.s,ylab="Score Values",main="Multivariate Score")
  abline(h=score1.s[cpm_s$changePoint],lty=2) 
}
return(out.data)
}


