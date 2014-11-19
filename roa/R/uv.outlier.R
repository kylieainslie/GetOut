#####################################
### Single Group Outlier Analysis ###
#####################################
uv.outlier<-function(data, annotation=NULL,annID, annName,cpmtype="Mann-Whitney",
                      save.path=NULL,diff=T,num.id=T,cut=0.85,common.genes=F,under=F,p=NULL){

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
num = dim(data)[1]

### GTI
  gti = function(data) (sum(as.numeric(data)>as.numeric((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T)))/num)*(mean(subset(data,data>((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T))))+((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T)))/mean(subset(data,data>((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T))))
  gti_under = function(data) (sum(as.numeric(data)<as.numeric((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T)))/num)*(mean(subset(data,data<((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T))))-((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T)))/mean(subset(data,data<((quantile(data,0.25,na.rm=T)-quantile(data,0.75,na.rm=T))+quantile(data,0.25,na.rm=T))))
### COPA
  copa = function(data) (quantile(data,cut,na.rm=T)-median(data))/mad(data)
  
### OSS
  oss_new = function(data) sum(subset(data,data>((quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))+quantile(data,0.75,na.rm=T)))-median(data))/mad(data)
  oss_under = function(data) sum(subset(data,data<(quantile(data,0.25,na.rm=T)-(quantile(data,0.75,na.rm=T)-quantile(data,0.25,na.rm=T))))-median(data))/mad(data)
### Variability
  variability = function(data) var(data)
  
#Convert data frame into a matrix
  data1<-data.matrix(data)
  
#Apply statistics
### GTI ############################################################################################
  if(under==F){data2 <- apply_pb(data1,2,gti,Title="Calculating GTI ...")}
  if(under==T){data2 <- apply_pb(data1,2,gti_under,Title="Calculating GTI ...")}
  data2.a<-subset(data2,is.nan(data2)==F)
  data2.b<-subset(data2.a,is.na(data2.a)==F)
  data2.c<-subset(data2.b,data2.b!="Inf")

if(is.null(cpmtype) & !is.null(p)){
  if(under==F){data3<-ifelse(data2.c>quantile(data2.c, p,na.rm=TRUE), data2.c, NA)}
  if(under==T){data3<-ifelse(data2.c<quantile(data2.c, p,na.rm=TRUE), data2.c, NA)}
}
if(!is.null(cpmtype) & is.null(p)){
  if(under==F){data2.s<-sort(data2,decreasing=T)}
  if(under==T){data2.s<-sort(data2)}
  
  if (diff==F){
    cpm_g<-detectChangePoint(data2.s,cpmType=cpmtype, ARL0=500, startup=20)
  }
  if (diff==T){
  diffg<-c(rep(0,length(data2.s)))
  for (i in 1:length(data2.s)-1){
    diffg[i]<-data2.s[i]-data2.s[i+1]
  }
  diffg1<-subset(diffg,is.nan(diffg)==F)
  diffg2<-subset(diffg1,is.na(diffg1)==F)
  diffg3<-subset(diffg2,diffg2!="Inf")
  cpm_g<-detectChangePoint(diffg3,cpmType=cpmtype, ARL0=500, startup=20)
  }
  if(under==F){data3<-ifelse(data2>=data2.s[cpm_g$changePoint], data2, NA)}
  if(under==T){data3<-ifelse(data2<=data2.s[cpm_g$changePoint], data2, NA)}
}  
  data4<-subset(data3, data3!="NA")
  
  if(num.id==T){
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
  if(num.id==F){
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
  write.csv(gti_out,file=paste(save.path,"gti_outlier.csv"))
  }
  #if (outsample==T & !is.null(save.path)){
  #  if(under==F){outlying.samples(data,names(data4),over_cut=cut,stat="GTI",savepath=save.path)}
  #  if(under==T){outlying.samples(data,names(data4),under_cut=cut,stat="GTI",savepath=save.path)}
  #}

  
#### COPA ##########################################################################################
  data2c <- apply_pb(data1,2,copa,Title="Calculating COPA ...")
  data2c.a<-subset(data2c,is.nan(data2c)==F)
  data2c.b<-subset(data2c.a,is.na(data2c.a)==F)
  data2c.c<-subset(data2c.b,data2c.b!="Inf")

if(is.null(cpmtype) & !is.null(p)){
  if(under==F){data3c<-ifelse(data2c.c>quantile(data2c.c, p,na.rm=TRUE), data2c.c, NA)}
  if(under==T){data3c<-ifelse(data2c.c<quantile(data2c.c, p,na.rm=TRUE), data2c.c, NA)}
}
if(!is.null(cpmtype) & is.null(p)){
  if(under==F){data2c.s<-sort(data2c,decreasing=T)}
  if(under==T){data2c.s<-sort(data2c)}

  if (diff==F){
    cpm_c<-detectChangePoint(data2c.s,cpmType=cpmtype, ARL0=500, startup=20)
  }
  if(diff==T){
  diffc<-c(rep(0,length(data2c.s)))
  for (i in 1:length(data2c.s)-1){
    diffc[i]<-data2c.s[i]-data2c.s[i+1]
  }
  diffc1<-subset(diffc,is.nan(diffc)==F)
  diffc2<-subset(diffc1,is.na(diffc1)==F)
  diffc3<-subset(diffc2,diffc2!="Inf")
  cpm_c<-detectChangePoint(diffc3,cpmType=cpmtype, ARL0=500, startup=20)
  }
  if(under==F){data3c<-ifelse(data2c>=data2c.s[cpm_c$changePoint], data2c, NA)}
  if(under==T){data3c<-ifelse(data2c<=data2c.s[cpm_c$changePoint], data2c, NA)}
}
data4c<-subset(data3c, data3c!="NA")
  
  if(num.id==T){
    noX<-substr(names(data4c), 2,15)
    noX1<-unique(noX)
    if(!is.null(annotation)){
    gsc<-annotation[annotation[,annID] %in% noX1,c(annID,annName)]
    Val<-data4c[noX1 %in% gsc[,1]]
    val<-Val[unique(names(Val))]
    }
    else {gsc<-noX1
          val<-as.vector(data4c)}
  }
  if(num.id==F){
    if(!is.null(annotation)){
    gsc<-annotation[annotation[,annID] %in% names(data4c),c(annID,annName)]
    val<-data4c[names(data4c) %in% gsc[,1]]  
    }
    else {gsc<-names(data4c)
          val<-as.vector(data4c)}
  }
  copa_out<-data.frame(gsc,Value=val,Stat="COPA")
  rownames(copa_out)<-NULL
  
  if(!is.null(save.path)){
  write.csv(copa_out,file=paste(save.path,"copa_outlier.csv"))
  }
  #if (outsample==T & !is.null(save.path)){
  #  if(under==F){outlying.samples(data,names(data4c),over_cut=cut,stat="copa",savepath=save.path)}
  #  if(under==T){outlying.samples(data,names(data4c),under_cut=cut,stat="copa",savepath=save.path)}  
  #}

  
### OSS #########################################################################################
  if(under==F){data2o <- apply_pb(data1,2,oss_new,Title="Calculating OSS ...")}
  if(under==T){data2o <- apply_pb(data1,2,oss_under,Title="Calculating OSS ...")}
    data2o.a<-subset(data2o,is.nan(data2o)==F)
    data2o.b<-subset(data2o.a,is.na(data2o.a)==F)
    data2o.c<-subset(data2o.b,data2o.b!="Inf")

if(is.null(cpmtype) & !is.null(p)){
  if(under==F){data3o<-ifelse(data2o.c>quantile(data2o.c, p,na.rm=TRUE), data2o.c, NA)}
  if(under==T){data3o<-ifelse(data2o.c<quantile(data2o.c, p,na.rm=TRUE), data2o.c, NA)}
}
if(!is.null(cpmtype) & is.null(p)){
  if(under==F){data2o.s<-sort(data2o.c,decreasing=T)}
  if(under==T){data2o.s<-sort(data2o.c)}

  if (diff==F){
    cpm_o<-detectChangePoint(data2o.s,cpmType=cpmtype, ARL0=500, startup=20)
  }
  
  if(diff==T){
  diffo<-c(rep(0,length(data2o.s)))
    for (i in 1:length(data2o.s)-1){
    diffo[i]<-data2o.s[i]-data2o.s[i+1]
    }
  diffo1<-subset(diffo,is.nan(diffo)==F)
  diffo2<-subset(diffo1,is.na(diffo1)==F)
  diffo3<-subset(diffo2,diffo2!="Inf")
  cpm_o<-detectChangePoint(diffo3,cpmType=cpmtype, ARL0=500, startup=20)
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
  write.csv(oss_out,file=paste(save.path,"oss_outlier.csv"))
  }
  #if (outsample==T & !is.null(save.path)){
  #  if(under==F){outlying.samples(data,names(data4o),over_cut=cut,stat="OSS",savepath=save.path)}
  #  if(under==T){outlying.samples(data,names(data4o),under_cut=cut,stat="OSS",savepath=save.path)}
  #}
    
  
### Variance #################################################################################
  data2v <- apply_pb(data1,2,variability,Title="Calculating variance ...")
  data2v.a<-subset(data2v,is.nan(data2v)==F)
  data2v.b<-subset(data2v.a,is.na(data2v.a)==F)
  data2v.c<-subset(data2v.b,data2v.b!="Inf")
  data2v.s<-sort(data2v.c,decreasing=T)
  
if(is.null(cpmtype) & !is.null(p)){
  data3v<-ifelse(data2v.c>quantile(data2v.c, p,na.rm=TRUE), data2v.c, NA)
}
if(!is.null(cpmtype) & is.null(p)){
  if (diff==F){
    cpm_v<-detectChangePoint(data2v.s,cpmType=cpmtype, ARL0=500, startup=20)
  }
  if (diff==T){
  diffv<-c(rep(0,length(data2v.s)))
  for (i in 1:length(data2v.s)-1){
    diffv[i]<-data2v.s[i]-data2v.s[i+1]
  }
  cpm_v<-detectChangePoint(diffv,cpmType=cpmtype, ARL0=500, startup=20)
  }
  data3v<-ifelse(data2v.c>=data2v.s[cpm_v$changePoint], data2v.c, NA)
}  
data4v<-subset(data3v, data3v!="NA")

  if(num.id==T){
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
  if(num.id==F){
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
  write.csv(var_out,file=paste(save.path,"var_outlier.csv"))
  }
  #if (outsample==T & !is.null(save.path)){
  #  outlying.samples(data,names(data4v),over_cut=cut,stat="var",savepath=save.path)
  #}
  
########################
### all results file ###
########################

if(!is.null(annotation)){
    all_results<-rbind(gti_out,copa_out,oss_out,var_out)
}
if(is.null(annotation)){all_results<-list(GTI=gti_out,COPA=copa_out,OSS=oss_out,Variance=var_out)}

if (!is.null(save.path)){
  write.csv(all_results,file=paste(save.path,"all_results.csv"))
}  
#################  
### cpm plots ###
#################
if(!is.null(save.path)){  
  setwd(save.path)
  par(mfrow=c(1,2))
### GTI
  jpeg(filename="gti_cpm_plots.jpg")
  cpm_plot_D_g<-plot(cpm_g$Ds)
  cpm_thresh_plot<-plot(sort(data2))
  abline(h=data2.s[cpm_g$changePoint],lty=2)
  dev.off()
### COPA
  jpeg(filename="copa_cpm_plots.jpg")
  cpm_plot_D_c<-plot(cpm_c$Ds)
  cpm_thresh_plot_c<-plot(sort(data2c))
  abline(h=data2c.s[cpm_c$changePoint],lty=2)
  dev.off()
### OSS
  jpeg(filename="oss_cpm_plots.jpg")
  cpm_plot_D_o<-plot(cpm_o$Ds)
  cpm_thresh_plot_o<-plot(sort(data2o))
  abline(h=data2o.s[cpm_o$changePoint],lty=2)
  dev.off()
### Var
  jpeg(filename="var_cpm_plots.jpg")
  cpm_plot_D_v<-plot(cpm_v$Ds)
  cpm_thresh_plot_v<-plot(sort(data2v))
  abline(h=data2v.s[cpm_v$changePoint],lty=2)
  dev.off()
}
if(is.null(save.path)){
  par(mfrow=c(2,2))
# GTI  
  cpm_thresh_plot<-plot(sort(data2),main="GTI",ylab="GTI Values")
  abline(h=data2.s[cpm_g$changePoint],lty=2)
# COPA  
  cpm_thresh_plot_c<-plot(sort(data2c),main="COPA",ylab="COPA Values")
  abline(h=data2c.s[cpm_c$changePoint],lty=2)
# OSS
  cpm_thresh_plot_o<-plot(sort(data2o),main="OSS",ylab="OSS Values")
  abline(h=data2o.s[cpm_o$changePoint],lty=2)
# Variance
  cpm_thresh_plot_v<-plot(sort(data2v),main="Variance",ylab="Variance values")
  abline(h=data2v.s[cpm_v$changePoint],lty=2)
}
par(mfrow=c(1,1))

#######################
### Genes in common ###
#######################
if(common.genes==T){
  g<-names(data4)
  c<-names(data4c)
  o<-names(data4o)
  v<-names(data4v)
  
  common4<-Reduce(intersect, list(g,c,o,v))
  if (length(common4)>0){
    method=5
  }
  if (length(common4)==0){
    warning("No outliers in common across all four methods")
    method=3
  }
  if (method==5){
    #Four methods
    if(num.id==T){common_genes<-substr(common4, 2,15)}
    else {common_genes<-common4}
    if(!is.null(annotation)){
    out4<-data.frame(annotation[annotation[,annID] %in% common_genes,],Methods=4)
    out4a<-out4[out4[,annName]!="---",]
    }
    else {out4a<-data.frame(common_genes,Methods=4)}
    
    #Three methods
    common3<-list(GCO=Reduce(intersect, list(g,c,o)),GCV=Reduce(intersect, list(g,c,v)),
                  GOV=Reduce(intersect, list(g,v,o)),COV=Reduce(intersect, list(v,c,o)))
    if(num.id==T){common_genes<-substr(unlist(common3), 2,15)}
    else {common_genes<-unlist(common3)}
    if(!is.null(annotation)){
    out3<-data.frame(annotation[annotation[,annID] %in% common_genes,],Methods=3)
    out3a<-out3[out3[,annName]!="---",]
    }
    else {out3a<-data.frame(common_genes, Methods=3)}
    common.results<-rbind(out4a,out3a)

  }
  if (method==3){
    common3<-list(GCO=Reduce(intersect, list(g,c,o)),GCV=Reduce(intersect, list(g,c,v)),
                 GOV=Reduce(intersect, list(g,v,o)),COV=Reduce(intersect, list(v,c,o)))
    common3a<-unique(as.vector(unlist(common3)))
    if (length(common3a)==0){
      common.genes=F
      warning("No outliers in common across three of the four statistics")}
    if(length(common3a)>0){
      if(num.id==T){common_genes<-substr(unlist(common3a), 2,15)}
      else {common_genes<-common3a}
        if(!is.null(annotation)){
        out3<-data.frame(annotation[annotation[,annID] %in% common_genes,],Methods=3)
        out3a<-out3[out3[,annName]!="---",]
        common.results<-out3a
        }
        else {out3a<-data.frame(common_genes, Methods=3)
              common.results<-out3a}
    
    }
  }  
}

##############
### Output ###
##############
if(common.genes==F){return(list(Outliers=all_results))}
if(common.genes==T){return(list(Outliers=all_results,Common_Genes=common.results))} 
}
