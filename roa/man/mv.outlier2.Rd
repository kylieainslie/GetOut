\name{mv.outlier2}
\alias{mv.outlier2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two Group Multivariate Outlier Statistic
}
\description{
Outliers can have biological relevance when determined using gene expression, RNA-Seq, methylation, or other types of data. Three approaches have been previously developed: Cancer Profile Outlier Analysis (COPA), Outlier Sum Statistic (OSS), and  Gene Tissue Index (GTI). This multivariate statistic combines COPA, OSS, and GTI into a score statistic to determine outliers.
}
\usage{
mv.outlier2(data_d, data_n, data_all, annotation = NULL, annID, annName, save.path = NULL, cut = 0.85, cpmtype = "Mann-Whitney", num.id = TRUE, under = FALSE, diff = TRUE, p = NULL)
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
  \item{cut}{
A number between 0 and 1 to be used as the cut-off for COPA calculation and the percentile above which outlying samples are determined.
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
  \item{diff}{
Logical indicating if consecutive differences should be calculated for the vector of ordered score statistics. If cpmtype is specified, the change point model would be applied to the vector of consecutive differences.
  }
  \item{p}{
A number between 0 and 1 to be used as the user defined cut-off above which outliers are defined.
}
}

\value{
A list with components
  \item{Outliers}{List of three-column data frames showing the outlier genes from each statistic.}
  \item{Common_Genes}{List of genes in common among all four or three of the four outlier statistics.}
}

\references{
Mpindi, J.P. et al. 2011. GTI: A novel algorithm for identifying outlier gene expression profiles from integrated microarray datasets. PLoS One 6(2): e17259.
Ross, G.J. Parametric and nonparametric sequential change detection in R: the cpm package. Journal of Statistical Software  (forthcoming).
Tibshirani, R. and T. Hastie. 2007. Outlier sums for differential gene expression analysis. Biostatistics 8(1): 2-8.
Tomlins, S.A. et al. 2005. Recurrent fusion of TMPRSS2 and ETS transcription factor genes in prostate cancer. Science 310: 644-648

}
\author{
Kylie Ainslie and Jeanne Kowalski
}

\seealso{
\code{\link{detectChangePoint}}, \code{\link{copa}}, \code{\link{gti}}, \code{\link{oss}}, \code{\link{mv.outlier}}
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
  mv.outlier2(dis.data,norm.data,all.data,annotation=ann,annID=1,annName=2,num.id=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{outlier}
\keyword{two group}% __ONLY ONE__ keyword per line
\keyword{multivariate}
