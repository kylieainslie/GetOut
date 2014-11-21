\name{uv.outlier}
\alias{uv.outlier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate Outlier Analysis
}
\description{
Outliers can have biological relevance when determined using gene expression, RNA-Seq, methylation, or other types of data. Three approaches have been previously developed: Cancer Profile Outlier Analysis (COPA), Outlier Sum Statistic (OSS), and  Gene Tissue Index (GTI). This function separately determines univariate outliers using the three previously mentioned outlier statistics in addition to a variance approach in a single group of samples. Results can be intersected to find genes identified in all four approaches or three of four approaches.
}
\usage{
uv.outlier(data, annotation = NULL, annID, annName, cpmtype = "Mann-Whitney", save.path = NULL, diff = TRUE, num.id = TRUE, cut = 0.85, common.genes = FALSE, under = FALSE, p = NULL)
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
  \item{cpmtype}{
The type of CPM which is used. For possible arguments see \code{\link{detectChangePoint}}
}
  \item{save.path}{
A character string indicating the path where the outputted results file is to be saved.
}
  \item{diff}{
Logical indicating if consecutive differences between sorted outlier statistics should be calculated before a change point is identified.
}
  \item{num.id}{
Logical indicating if data column names begin with a number.
}
  \item{cut}{
A number between 0 and 1 to be used as the cut-off for COPA calculation and the percentile above which outlying samples are determined.
}
  \item{common.genes}{
Logical indicating if results from the four outlier approaches should be intersected.
}
  \item{under}{
Logical indicating if low-valued outliers should be identified.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{detectChangePoint}}, \code{\link{outlying.samples}}, \code{\link{copa}}, \code{\link{gti}}, \code{\link{oss}}, \code{\link{variance}}, \code{\link{mv.outlier}}, \code{\link{uv.outlier2}}
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
  
## Example with annotation file ##
  uv.outlier(Data,annotation=ann,annID=1,annName=2,num.id=FALSE, common.genes=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{outlier}
\keyword{univariate}% __ONLY ONE__ keyword per line
