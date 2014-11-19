################################################################################
### Two Group Outlier Analysis ###
################################################################################
uv.outlier2<-function(data_d,data_n,data_all, annotation=NULL,annID, annName, diff=T,cpmtype="Mann-Whitney",
                       save.path=NULL, cut=0.85,num.id=T,common.genes=T,under=F,p=NULL){

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
  if(under==T & cut>0.5){warning("cut value doesn't match under status")}
  if(under==F & cut<0.5){warning("cut value doesn't match under status")}
  
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
##Define number of samples per probe
  
  num1 = dim(data_n)[1]
  num2 = dim(data_d)[1]

#GTI - a statistic is calculated for each group (1=normal; 2=disease)
gti1 = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T)))/num1)*(mean(subset(data,data>((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T))))+((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T)))/mean(subset(data,data>((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T))))  
gti2 = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T)))/num2)*(mean(subset(data,data>((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T))))+((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T)))/mean(subset(data,data>((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T))))  

gti_under1 = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T)))/num1)*(mean(subset(data,data<((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T))))-((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T)))/mean(subset(data,data<((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T))))
gti_under2 = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T)))/num2)*(mean(subset(data,data<((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T))))-((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T)))/mean(subset(data,data<((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T))))

gti<-function(data_d,data_n,data_all,Under=under){
  if(Under==F){
    data_norm <- apply_pb(data1,2,gti1,Title="Calculating GTI ...")
    data_dis<- apply_pb(data2,2,gti2,Title="Calculating GTI ...")
  }
  if(Under==T){
    data_norm <- apply_pb(data1,2,gti_under1,Title="Calculating GTI ...")
    data_dis<- apply_pb(data2,2,gti_under2,Title="Calculating GTI ...")
  }
  dist<-data_dis-data_norm
return(dist)  
}
#COPA             
# use combined (normal and disease) samples for median and mad

#take rth percentile of disease group
copa<-function(data_d,data_n,data_all,Cut=0.85){
  #calculate quantile specified by cut for each gene
    r<-function(data) (quantile(data,Cut,na.rm=T)) #quantile function
    x<-data.matrix(data_d)
    x1<-apply_pb(x,2,r)
  #calculate median
    #median1<- function(data) quantile(data,0.5,na.rm=T)
    data1<-data.matrix(data_all)
    med<-apply_pb(data1,2,median,Title="Calculating median ...")
  #calculate mad
    #mad1<-function(data) mad(data,na.rm=T)
    mad1<-apply_pb(data1,2,mad,Title="Calculating median absolute deviations ...")
  #calculate copa
    copa<-(x1-med)/mad1
return(copa)
}

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

dat.gti<-gti(data1,data2,data3)
dat.copa<-copa(data1,data2,data3,Cut=cut)
dat.oss<-oss(data1,data2,data3)

#####################################
### Change point model and output ###
#####################################

### GTI
if(is.null(cpmtype) & !is.null(p)){
  if(under==F){dat.gti3<-ifelse(dat.gti>quantile(dat.gti, p,na.rm=TRUE), dat.gti, NA)}
  if(under==T){dat.gti3<-ifelse(dat.gti<quantile(dat.gti, p,na.rm=TRUE), dat.gti, NA)}
}

if(!is.null(cpmtype) & is.null(p)){
  if(under==F){dat.gti.s<-sort(dat.gti,decreasing=T)}
  if(under==T){dat.gti.s<-sort(dat.gti)}

  #determining change point on ordered raw data (highest to lowest)
  if (diff==F){cpm_g<-detectChangePoint(dat.gti.s,cpmType=cpmtype, ARL0=500, startup=20)}
  
  #determine change point for difference
  if (diff==T){
    diffg<-c(rep(0,length(dat.gti)))
    for (i in 1:length(dat.gti)-1){
      diffg[i]<-dat.gti.s[i]-dat.gti.s[i+1]
    }
    diffg1<-subset(diffg,is.nan(diffg)==F)
    diffg2<-subset(diffg1,is.na(diffg1)==F)
    diffg3<-subset(diffg2,diffg2!="Inf")

    cpm_g<-detectChangePoint(diffg3,cpmType=cpmtype, ARL0=500, startup=20)
  }  

  dat.gti2<-ifelse(dat.gti>=dat.gti.s[cpm_g$changePoint], dat.gti, NA)
  dat.gti3<-subset(dat.gti2, dat.gti2!="NA")
}
dat.gti4<-dat.gti3[unique(names(dat.gti3))]
  
  if(num.id==T){
    noX<-substr(names(dat.gti4), 2,15)
    noX1<-unique(noX)
    if(!is.null(annotation)){
      gs<-annotation[annotation[,annID] %in% noX1,c(annID,annName)]
      Val<-dat.gti4[noX1 %in% gs[,1]]
      val<-Val[unique(names(Val))]
    }
    else {gs<-noX1
          val<-as.vector(dat.gti4)}
  }
  if(num.id==F){
    if(!is.null(annotation)){
      gs<-annotation[annotation[,annID] %in% names(dat.gti4),c(annID,annName)]
      val<-dat.gti4[names(dat.gti4) %in% gs[,1]]
    }
    else {gs<-names(dat.gti4)
          val<-as.vector(dat.gti4)} 
  }
  gti_out<-data.frame(gs, Value=val,Stat="GTI")
  rownames(gti_out)<-NULL
  write.csv(gti_out,file=paste(save.path,"gti_outlier2.csv"))

### COPA
if(is.null(cpmtype) & !is.null(p)){
  if(under==F){dat.copa3<-ifelse(dat.copa>quantile(dat.copa, p,na.rm=TRUE), dat.copa, NA)}
  if(under==T){dat.copa3<-ifelse(dat.copa<quantile(dat.copa, p,na.rm=TRUE), dat.copa, NA)}
}

if(!is.null(cpmtype) & is.null(p)){
  if(under==F){dat.copa.s<-sort(dat.copa,decreasing=T)}
  if(under==T){dat.copa.s<-sort(dat.copa)}

  #determining change point on ordered raw data (highest to lowest)
  if (diff==F){cpm_c<-detectChangePoint(dat.copa.s,cpmType=cpmtype, ARL0=500, startup=20)}

  #determine change point for difference
  if (diff==T){
    diffc<-c(rep(0,length(dat.copa)))
  for (i in 1:length(dat.copa)-1){
    diffc[i]<-dat.copa.s[i]-dat.copa.s[i+1]
  }
  diffc1<-subset(diffc,is.nan(diffc)==F)
  diffc2<-subset(diffc1,is.na(diffc1)==F)
  diffc3<-subset(diffc2,diffc2!="Inf")
  
  cpm_c<-detectChangePoint(diffc3,cpmType=cpmtype, ARL0=500, startup=20)
}  

  dat.copa2<-ifelse(dat.copa>=dat.copa.s[cpm_c$changePoint], dat.copa, NA)
  dat.copa3<-subset(dat.copa2, dat.copa2!="NA")
}
dat.copa4<-dat.copa3[unique(names(dat.copa3))]

if(num.id==T){
  noX<-substr(names(dat.copa4), 2,15)
  noX1<-unique(noX)
  if(!is.null(annotation)){
    gs<-annotation[annotation[,annID] %in% noX1,c(annID,annName)]
    Val<-dat.copa4[noX1 %in% gs[,1]]
    val<-Val[unique(names(Val))]
  }
  else {gs<-noX1
        val<-as.vector(dat.copa4)}
}
if(num.id==F){
  if(!is.null(annotation)){
    gs<-annotation[annotation[,annID] %in% names(dat.copa4),c(annID,annName)]
    val<-dat.copa4[names(dat.copa4) %in% gs[,1]]
  }
  else {gs<-names(dat.copa4)
        val<-as.vector(dat.copa4)}
}
copa_out<-data.frame(gs, Value=val,Stat="COPA")
rownames(copa_out)<-NULL
write.csv(copa_out,file=paste(save.path,"copa_outlier2.csv"))


### OSS
if(is.null(cpmtype) & !is.null(p)){
  if(under==F){dat.oss3<-ifelse(dat.oss>quantile(dat.oss, p,na.rm=TRUE), dat.oss, NA)}
  if(under==T){dat.oss3<-ifelse(dat.oss<quantile(dat.oss, p,na.rm=TRUE), dat.oss, NA)}
}

if(!is.null(cpmtype) & is.null(p)){
  if(under==F){dat.oss.s<-sort(dat.oss,decreasing=T)}
  if(under==T){dat.oss.s<-sort(dat.oss)}

  #determining change point on ordered raw data (highest to lowest)
  if (diff==F){cpm_o<-detectChangePoint(dat.oss.s,cpmType=cpmtype, ARL0=500, startup=20)}

  #determine change point for difference
  if (diff==T){
  diffo<-c(rep(0,length(dat.oss)))
  for (i in 1:length(dat.oss)-1){
    diffo[i]<-dat.oss.s[i]-dat.oss.s[i+1]
  }
  diffo1<-subset(diffo,is.nan(diffo)==F)
  diffo2<-subset(diffo1,is.na(diffo1)==F)
  diffo3<-subset(diffo2,diffo2!="Inf")
  
  cpm_o<-detectChangePoint(diffo3,cpmType=cpmtype, ARL0=500, startup=20)
}  

  dat.oss2<-ifelse(dat.oss>=dat.oss.s[cpm_o$changePoint], dat.oss, NA)
  dat.oss3<-subset(dat.oss2, dat.oss2!="NA")
}
dat.oss4<-dat.oss3[unique(names(dat.oss3))]

if(num.id==T){
  noX<-substr(names(dat.oss4), 2,15)
  noX1<-unique(noX)
  if(!is.null(annotation)){
    gs<-annotation[annotation[,annID] %in% noX1,c(annID,annName)]
    Val<-dat.oss4[noX1 %in% gs[,1]]
    val<-Val[unique(names(Val))]
  }
  else {gs<-noX1
        val<-as.vector(dat.oss4)}
}
if(num.id==F){
  if(!is.null(annotation)){
    gs<-annotation[annotation[,annID] %in% names(dat.oss4),c(annID,annName)]
    val<-dat.oss4[names(dat.oss4) %in% gs[,1]]
  }
  else {gs<-names(dat.oss4)
        val<-as.vector(dat.oss4)}
}
oss_out<-data.frame(gs, Value=val,Stat="OSS")
rownames(oss_out)<-NULL
write.csv(oss_out,file=paste(save.path,"oss_outlier2.csv"))
#if (outsample==T){
#  if(under==F){outlying.samples(data_all,names(dat.oss4),over_cut=cut,stat="OSS",savepath=save.path)}
#  if(under==T){outlying.samples(data_all,names(dat.oss4),under_cut=cut,stat="OSS",savepath=save.path)}
#}
##############
### Output ###
##############
if(!is.null(annotation)){
  all_results<-rbind(gti_out,copa_out,oss_out)
}
if(is.null(annotation)){all_results<-list(GTI=gti_out,COPA=copa_out,OSS=oss_out)}

if (!is.null(save.path)){
  write.csv(all_results,file=paste(save.path,"all_results.csv"))
}
#################
### cpm plots ###
#################
if(!is.null(save.path)){
setwd(save.path)
# GTI
jpeg(filename="out_plot_gti.jpg")

par(mfrow=c(1,2))
cpm_plot_D_g<-plot(cpm_g$Ds)
cpm_thresh_plot_g<-plot(sort(dat.gti))
abline(h=dat.gti.s[cpm_g$changePoint],lty=2)

dev.off()
# COPA
jpeg(filename="out_plot_copa.jpg")

par(mfrow=c(1,2))
cpm_plot_D_c<-plot(cpm_c$Ds)
cpm_thresh_plot_c<-plot(sort(dat.copa))
abline(h=dat.copa.s[cpm_c$changePoint],lty=2)

dev.off()
# OSS
jpeg(filename="out_plot_oss.jpg")

par(mfrow=c(1,2))
cpm_plot_D_o<-plot(cpm_o$Ds)
cpm_thresh_plot_o<-plot(sort(dat.oss))
abline(h=dat.oss.s[cpm_o$changePoint],lty=2)

dev.off()
#reset window
par(mfrow=c(1,1))
}
if(is.null(save.path)){
  par(mfrow=c(3,1))
# GTI  
  cpm_thresh_plot_g<-plot(sort(dat.gti),main="GTI",ylab="GTI Values")
  abline(h=dat.gti.s[cpm_g$changePoint],lty=2)
# COPA
  cpm_thresh_plot_c<-plot(sort(dat.copa),main="COPA",ylab="COPA Values")
  abline(h=dat.copa.s[cpm_c$changePoint],lty=2)
# OSS
  cpm_thresh_plot_o<-plot(sort(dat.oss),main="OSS",ylab="OSS Values")
  abline(h=dat.oss.s[cpm_o$changePoint],lty=2)

par(mfrow=c(1,1))
}
#######################
### Genes in common ###
#######################
if(common.genes==T){
  g<-names(dat.gti4)
  c<-names(dat.copa4)
  o<-names(dat.oss4)
  
  common3<-Reduce(intersect, list(g,c,o))

  if (length(common3)==0){
    warning("No outliers in common across all three methods")
    common.genes=F
  }
  if (length(common3)>0){
    #Four methods
    if(num.id==T){common_genes<-substr(common3, 2,15)}
    else {common_genes<-common3}
    if(!is.null(annotation)){
      out3<-data.frame(annotation[annotation[,annID] %in% common_genes,],Methods=3)
      out3a<-out3[out3[,annName]!="---",]
    }
    else {out3a<-data.frame(common_genes,Methods=3)}

    common.results<-out3a
    
  }  
}

##############
### Output ###
##############
if(common.genes==F){return(list(Outliers=all_results))}
if(common.genes==T){return(list(Outliers=all_results,Common_Genes=common.results))} 
}

 
