\name{mv.outlier}
\alias{mv.outlier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Outlier Statistic
}
\description{
Outliers can have biological relevance when determined using gene expression, RNA-Seq, methylation, or other types of data. Three approaches have been previously developed: Cancer Profile Outlier Analysis (COPA), Outlier Sum Statistic (OSS), and  Gene Tissue Index (GTI). This multivariate statistic combines COPA, OSS, and GTI into a score statistic to determine outliers.
}
\usage{
mv.outlier(data, annotation=NULL, annID, annName, save.path=NULL, cut=0.85,cpmtype="Mann-Whitney",num.id=T, under=F, diff=T, p=NULL)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{cut}{
A number between 0 and 1 that indicates the percentile above which samples are determined to be "outlying" when calculating the COPA statistic.
}
  \item{cpmtype}{
The type of CPM which is used. For possible arguments see \code{\link{detectChangePoint}}
}
  \item{num.id}{
Logical indicating if the IDs (and column names of data) start with a number. R adds an "X" to column names that begin with a number, so the "X" is removed in the results file if num.id=T.
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
A three-column data frame providing the gene symbol and value of the score statistic for each ID considered an outlier.
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
\code{\link{detectChangePoint}}, \code{\link{outlying.samples}}, \code{\link{copa}}, \code{\link{gti}}, \code{\link{oss}}
}
\examples{
#generate non-outlier data
 sim.data<-matrix(c(rnorm(1000*200,0,1)),nrow=200)
  
#generate outlier gene expression
  y<-matrix(rnorm(200*10,0,1),nrow=200)
  m<-max(sim.data,y)
  c<-c(rep(m,20),rep(0,180))
  test.genes<-y+c
  Data<-cbind(test.genes,sim.data)
  
#add column names
  pre<-"test.gene"
  suf<-seq(1:10)
  prefix<-"gene"
  suffix<-seq(1:1000)
  colnames(Data)<-c(paste(pre,suf,sep=""),paste(prefix,suffix,sep=""))

#create annotation file
  p<-"gene.symbol"
  s<-seq(1:1010)
  ann<-data.frame(ID=colnames(Data),Gene=paste(p,s,sep=""))
  
  mv.outlier(Data,annotation=ann,annID=1,annName=2,num.id=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{outlier}
\keyword{multivariate}
