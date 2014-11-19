#########################################################
### Determining samples with outlying gene expression ###
#########################################################
outlying.samples<-function(data,outliers,cut=0.85,cpm=T,cpmtype="Mann-Whitney",
                           savepath=NULL,diff=T,under=F){  

if(cpm==T){cut=NULL}  
#extract probes from entire data set
out_probes<-data[,colnames(data) %in% outliers]
out_probes<-out_probes[,unique(colnames(out_probes))]

#define parameters                       
#cutoff<-c(rep(0,length(outliers))) #outlier cut-off
x<-dim(data)[1] #number of samples (rows)
y<-dim(out_probes)[2] #number of columns in outlier data set (number of outlier probes)

out_ind<-matrix(c(rep(0,length(outliers)*x)), nrow=x,ncol=length(outliers)) #indicator matrix

#determine which samples had outlying expression values
if(under==F){
  cutoff<-c(rep(0,y))
  if(!is.null(cut)){
    for (i in 1:y){
      cutoff[i]<-quantile(out_probes[,i],cut,na.rm=T)
      out_ind[,i]<-ifelse(out_probes[,i]>cutoff[i],1,0)
    }
  }
  if(cpm==T){
    for(i in 1:y){
      sort.out<-sort(out_probes[,i],decreasing=T)
      if (diff==F){cpm_s<-detectChangePoint(sort.out,cpmType=cpmtype,ARL0=500,startup=20)}
      if (diff==T){
        diffs<-c(rep(0,length(sort.out)))
        for (j in 1:length(sort.out)-1){
          diffs[j]<-sort.out[j]-sort.out[j+1]
        }
        #diffs1<-subset(diffs,is.nan(diffs)==F)
        #diffs2<-subset(diffs1,is.na(diffs1)==F)
        #diffs3<-subset(diffs2,diffs2!="Inf")
        
        cpm_s<-detectChangePoint(diffs,cpmType=cpmtype,ARL0=500,startup=20)
      }
      cutoff[i]<-sort.out[cpm_s$changePoint]
      out_ind[,i]<-ifelse(out_probes[,i]>cutoff[i],1,0)
    }
  } 
}
if(under==T){
  cutoff<-c(rep(0,y))
  if(!is.null(cut)){
    for (i in 1:y){
      cutoff[i]<-quantile(out_probes[,i],cut,na.rm=T)
      out_ind[,i]<-ifelse(out_probes[,i]<cutoff[i],1,0)
    }
  }
 if(cpm==T){
   for(i in 1:y){
     sort.out<-sort(out_probes[,i])
     cpm_s<-detectChangePoint(sort.out,cpmType=cpmtype, ARL0=500, startup=20)
     if (diff==F){cpm_s<-detectChangePoint(sort.out,cpmType=cpmtype, ARL0=500, startup=20)}
     if (diff==T){
       diffs<-c(rep(0,length(sort.out)))
       for (j in 1:length(sort.out)-1){
         diffs[j]<-sort.out[j]-sort.out[j+1]
       }
       diffs1<-subset(diffs,is.nan(diffs)==F)
       diffs2<-subset(diffs1,is.na(diffs1)==F)
       diffs3<-subset(diffs2,diffs2!="Inf")
       
       cpm_s<-detectChangePoint(diffs,cpmType=cpmtype, ARL0=500, startup=20)
     }
     cutoff[i]<-sort.out[cpm_s$changePoint]
     out_ind[,i]<-ifelse(out_probes[,i]<cutoff[i],1,0)
   }
 } 
}
#############
### Plots ###
#############
if(!is.null(savepath)){
  setwd(savepath)
  pdf(file="outlying_sample_plots2.pdf")
  
  for(i in 1:y){
    if(under==F){sort.out<-sort(out_probes[,i],decreasing=T)}
    if(under==T){sort.out<-sort(out_probes[,i])}
    plot(sort.out,ylab="Score Values",xlab="Sample",main=colnames(out_probes)[i])
    abline(h=cutoff[i],lty=2)
  } 
  dev.off()
}

#cbind sample IDs
sample_id<-data.frame(rownames(data),out_ind)
if(mode(outliers)=="numeric"){outliers<-as.character(outliers)}
  names(sample_id)<-c("ID",outliers)
 sample_ind<-c(rep(0,dim(sample_id)[1]))
 for(i in 1:dim(sample_id)[1]){
sample_ind[i]<-ifelse(sum(sample_id[i,-1])>0,1,0)
}
out_sample<-data.frame(rownames(data),sample_ind)
names(out_sample)<-c("ID","Outlying")

if(!is.null(savepath)){
  setwd(savepath)  
    write.csv(out_sample,file="outlying_sample_indicator.csv")
    write.csv(sample_id,file="outlying_samples_by_probe.csv")
  return("Results files have been output to specified path.")
}
if(is.null(savepath)){
  return(list(by.outlier=sample_id,sample.ind=out_sample))
}

}


       