##################################################
### Multivariate Outlier statistic - Two Group ###
##################################################
mv.outlier2<-function(data_d,data_n,data_all,annotation=NULL,annID,annName,save.path=NULL,cut=0.85,
                     cpmtype="Mann-Whitney",num.id=TRUE, under=FALSE, diff=TRUE, p=NULL){

  #check for missing arguments
  if(missing(data_d) | missing(data_n) | missing(data_all)){
    stop("No data specified")
  }
  #check for numeric arguments
  if(!is.numeric(data_d) | !is.numeric(data_n) | !is.numeric(data_all)){
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
#Define number of samples per probe
  
  num1 = dim(data_n)[1]
  num2 = dim(data_d)[1]
  
#GTI - a statistic is calculated for each group (1=normal; 2=disease)
  gti1 = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/num1)*(mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))+((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))  
  gti2 = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/num2)*(mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))+((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE)))/mean(subset(data,data>((quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))+quantile(data,0.75,na.rm=TRUE))))  

  gti_under1 = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/num1)*(mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))-((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))
  gti_under2 = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/num2)*(mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))-((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE)))/mean(subset(data,data<((quantile(data,0.25,na.rm=TRUE)-quantile(data,0.75,na.rm=TRUE))+quantile(data,0.25,na.rm=TRUE))))

  gti<-function(data_d,data_n,data_all,Under=under){
    if(Under==FALSE){
      data_norm <- apply_pb(data1,2,gti1,Title="Calculating GTI ...")
      data_dis<- apply_pb(data2,2,gti2,Title="Calculating GTI ...")
    }
    if(Under==TRUE){
      data_norm <- apply_pb(data1,2,gti_under1,Title="Calculating GTI ...")
      data_dis<- apply_pb(data2,2,gti_under2,Title="Calculating GTI ...")
    }
    dist<-data_dis-data_norm
    return(dist)  
  }
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
    med<-apply_pb(data1,2,median,Title="Calculating medians for COPA ...")
    #calculate mad
    #mad1<-function(data) mad(data,na.rm=TRUE)
    mad1<-apply_pb(data1,2,mad,Title="Calculating median absolute deviations for COPA ...")
    #calculate copa
    copa<-(x1-med)/mad1
    return(copa)
  }
  
#OSS
  oss<-function(data_d,data_n,data_all,Under=under){
    #calculate median
    data1<-data.matrix(data_all)
    med<-apply_pb(data1,2,median,Title="Calculating median for OSS ...")
    #calculate mad
    data2<-data.matrix(data_n)
    mad2<-apply_pb(data2,2,mad,Title="Calculating median absolute deviations for OSS ...")
    #standardize gene expression
    x2<-sweep(data_all,2,med,FUN="-")
    x3<-sweep(x2,2,mad2,FUN="/")
    #determine cutoff
    if(Under==FALSE){cutoff<-function(data) quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE)+quantile(data,0.75,na.rm=TRUE)}
    if(Under==TRUE){cutoff<-function(data) quantile(data,0.25,na.rm=TRUE)-(quantile(data,0.75,na.rm=TRUE)-quantile(data,0.25,na.rm=TRUE))}
    
    gene_cut<-apply_pb(data.matrix(x3),2,cutoff,Title="Calculating cut-offs for OSS ...")
    
    #standardize disease group gene expression
    y2<-sweep(data_d,2,med,FUN="-")
    y3<-sweep(x2,2,mad2,FUN="/")
    y4<-data.matrix(y3)   
    #calculate oss
    oss<-c(rep(0,dim(y4)[2]))
    for (i in 1:dim(y4)[2]){
      if(Under==FALSE){oss[i]<-sum(subset(y4[,i],y4[,i]>gene_cut[i]))}
      if(Under==TRUE){oss[i]<-sum(subset(y4[,i],y4[,i]<gene_cut[i]))}
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
  
  dat.gti<-gti(data1,data2,data3)
  dat.copa<-copa(data1,data2,data3,Cut=cut)
  dat.oss<-oss(data1,data2,data3)
  
### calculate score statistic
  score<-c(rep(0,length(dat.gti)))
  for(i in 1:length(dat.gti)){
    #print(i)
    value<-matrix(c(dat.gti[i],dat.copa[i],dat.oss[i]),nrow=3)
    score[i]<-t(value)%*%diag(3)%*%value
  }  
  score1<-data.matrix(score)
  names(score1)<-names(dat.gti)
  
  
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
      
      cpm_s<-detectChangePoint(diffs,cpmType=cpmtype, ARL0=500, startup=20)
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
    write.csv(out.data,paste("mv_outlier2_",cut,".csv",sep=""))
  }
#############
### Plots ###
#############
if(!is.null(save.path)){
  setwd(save.path)
  jpeg(filename=paste("mv_outlier2_plot",cut,".jpg",sep=""))
  
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

