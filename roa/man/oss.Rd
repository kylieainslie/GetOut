\name{oss}
\alias{oss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Outlier Sum Statistic
}
\description{
Outliers can have biological relevance when determined using gene expression, RNA-Seq, methylation, or other types of data. Tibshirani and Hastie introduced the Outlier Sum Statistic (OSS) as a method of detecting outliers in cancer samples. This function calculates OSS and determines outliers using a change point model.
}
\usage{
oss(data, annotation = NULL, annID, annName, save.path = NULL, diff = TRUE, cpmtype = "Mann-Whitney", num.id = TRUE, under = FALSE, p = NULL)
}
\arguments{
  \item{data}{
A numeric matrix of values arranged so that genes/IDs are columns and samples are rows. Gene/ID names should be column names.
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
The type of CPM which is used. For possible arguments see \code{\link{detectChangePoint}}
}
  \item{num.id}{
Logical indicating if data column names begin with a number.
}
  \item{under}{
Logical indicating whether low-valued outliers should be identified.
}
  \item{p}{
A number between 0 and 1 to be used as the user defined cut-off above which outliers are defined.
}
}

\value{
  A three- or four-column data frame showing the outlier genes identified by the OSS statistic.
}
\references{
Ross, G.J. Parametric and nonparametric sequential change detection in R: the cpm package. Journal of Statistical Software  (forthcoming).
Tibshirani, R. and T. Hastie. 2007. Outlier sums for differential gene expression analysis. Biostatistics 8(1): 2-8.
}
\author{
Kylie Ainslie and Jeanne Kowalski
}

\seealso{
\code{\link{detectChangePoint}}, \code{\link{uv.outlier}}, \code{\link{uv.outlier2}}, \code{\link{mv.outlier}}
}
\examples{
#generate non-outlier data
 sim.data<-matrix(c(rnorm(1000*200,0,1)),nrow=200)
  
#generate outlier gene expression
  y<-matrix(rnorm(200*10,0,1),nrow=200)
  m<-max(sim.data,y)
  c<-c(rep(m,20),rep(0,180))
  test.genes<-y+c
  under.genes<-y-c
  Data<-cbind(test.genes,sim.data)
  Data.under<-cbind(under.genes,sim.data)
  
#add column names
  pre<-"test.gene"
  pre.u<-"under.gene"
  suf<-seq(1:10)
  prefix<-"gene"
  suffix<-seq(1:1000)
  colnames(Data)<-c(paste(pre,suf,sep=""),paste(prefix,suffix,sep=""))
  colnames(Data.under)<-c(paste(pre.u,suf,sep=""),paste(prefix,suffix,sep=""))

## Example without annotation file ##
  oss(Data,num.id=FALSE) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{oss}
\keyword{outlier}% __ONLY ONE__ keyword per line
