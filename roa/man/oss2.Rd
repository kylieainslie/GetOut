\name{oss2}
\alias{oss2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two Group Outlier Sum Statistic
}
\description{
Outliers can have biological relevance when determined using gene expression, RNA-Seq, methylation, or other types of data. Tibshirani and Hastie introduced the Outlier Sum Statistic (OSS) as a method of detecting outliers in cancer samples versus normal samples. This function calculates the two group OSS and determines outliers using either a change point model or a user-defined threshold.
}
\usage{
oss2(data_d, data_n, data_all, annotation = NULL, annID, annName, save.path = NULL, diff = TRUE, cpmtype = "Mann-Whitney", num.id = TRUE, under = FALSE, p = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_d}{
A numeric matrix of values from samples in the disease group arranged so that genes/IDs are columns and samples are rows. Gene/ID names should be column names.
}
  \item{data_n}{
A numeric matrix of values from samples in the control group arranged so that genes/IDs are columns and samples are rows. Gene/ID names should be column names.
}
  \item{data_all}{
A numeric matrix of values from samples in both the disease group and control group arranged so that genes/IDs are columns and samples are rows. Gene/ID names should be column names.
}
  \item{annotation}{
A dataframe or matrix that contains a column of IDs and a column of the corresponding gene symbol for each ID.
}
  \item{annID}{
Positive integer that corresponds to the column in annotation that contains the IDs.
}
  \item{annName}{
Positive integer that corresponds to the column in annotation that contains the gene symbols.
}
  \item{save.path}{
A character string indicating the path where the outputted results file is to be saved.
}
  \item{diff}{
Logical indicating if consecutive differences between sorted outlier statistics should be calculated before a change point is identified.
}
  \item{cpmtype}{
The type of CPM which is used. For possible arguments see \code{\link{detectChangePoint}}.
}
  \item{num.id}{
Logical indicating if data column names begin with a number.
}
  \item{under}{
Logical indicating if low-valued outliers should be identified.
  }
  \item{p}{
A number between 0 and 1 to be used as the user defined cut-off above which outliers are defined.
}
}

\value{
A three- or four-column data frame showing the outlier genes identified by the OS statistic.
}
\references{
Ross, G.J. Parametric and nonparametric sequential change detection in R: the cpm package. Journal of Statistical Software  (forthcoming).
Tibshirani, R. and T. Hastie. 2007. Outlier sums for differential gene expression analysis. Biostatistics 8(1): 2-8.
}
\author{
Kylie Ainslie and Jeanne Kowalski
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{detectChangePoint}}, \code{\link{uv.outlier2}}, \code{\link{mv.outlier2}}
}
\examples{
#generate normal data
 norm.data<-matrix(c(rnorm(1010*100,0,1)),nrow=100)
  
#generate disease data with 10 outlier genes
  mydata<-matrix(rnorm(1000*100,0,1),nrow=100)
  y<-matrix(rnorm(100*10,0,1),nrow=100)
  m<-max(mydata,y)
  c<-c(rep(m,20),rep(0,80))
  test.genes<-y+c
  dis.data<-cbind(test.genes,mydata)

#add column names
  pre<-"test.gene"
  suf<-seq(1:10)
  prefix<-"gene"
  suffix<-seq(1:1000)
  colnames(dis.data)<-c(paste(pre,suf,sep=""),paste(prefix,suffix,sep=""))
  colnames(norm.data)<-c(paste(pre,suf,sep=""),paste(prefix,suffix,sep=""))
  
#complete data
  all.data<-rbind(dis.data,norm.data)
  rownames(all.data)<-paste("sample",seq(1:200),sep="")

#create annotation file
  p<-"gene.symbol"
  s<-seq(1:1010)
  ann<-data.frame(ID=colnames(dis.data),Gene=paste(p,s,sep=""))
  
## Example with annotation file ##
  oss2(dis.data,norm.data,all.data,annotation=ann,annID=1,annName=2,num.id=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{oss}
\keyword{outlier}% __ONLY ONE__ keyword per line
