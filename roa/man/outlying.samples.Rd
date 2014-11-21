\name{outlying.samples}
\alias{outlying.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determining Outlying Samples
}
\description{
When outliers are identified among samples, it is often of interest to find out which samples have outlying expression for each outlier. This function identifies outlying samples for each outlier and samples that have outlying values for at least one outlier.
}
\usage{
outlying.samples(data, outliers, cut = 0.85,cpm = TRUE, cpmtype = "Mann-Whitney", savepath = NULL, diff = TRUE, under = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A numeric matrix of values arranged so that genes/IDs are columns and samples are rows. Gene/ID names should be column names.
}
  \item{outliers}{
A vector of column names that identify the outlyier columns from data.
}
  \item{cut}{
A number between 0 and 1 indicating the cut-off by which outlying sample values are determined.
}
  \item{cpm}{
Logical indicating if a change point model should be used to identify outlying samples.
}
  \item{cpmtype}{
The type of change point model to be used to identify outlying samples. For possible arguments see \code{\link{detectChangePoint}}
}
  \item{savepath}{
A character string indicating the path where the outputted results file is to be saved. If save.path is specified, results are outputted in a .csv file.
}
  \item{diff}{
Logical indicating if consecutive differences should be calculated for the vector of ordered score statistics. If cpmtype is specified, the change point model would be applied to the vector of consecutive differences.
}
  \item{under}{
Logical indicating if outliers are to be identifed above over_cut or below under_cut.
}
}
\value{
If savepath is not NULL, two separate .csv files will be created for each statistic: 1) an indicator for each outlier indicating those samples that have outlying values and 2) an indicator of which samples have outlying values for at least one outlier. If savepath=NULL a list with the following components is output:
  \item{by.outlier}{A matrix with each row representing a sample and each column representing an outlier. A 1 in element [i,j] indicates that sample i has an outlying value (as detrmined by a change point model or cut) for outlier j. }
  \item{sample.ind}{A two column matrix with the sample IDs in column one and an indicator of whether each sample has an outlying value in at least one outlier.}

}

\author{
Kylie Ainslie and Jeanne Kowalski
}

\examples{
#generate non-outlier data
 sim.data<-matrix(c(rnorm(1000*100,0,1)),nrow=100)
  
#generate outlier gene expression
  y<-matrix(rnorm(100*10,0,1),nrow=100)
  m<-max(sim.data,y)
  c<-c(rep(m,10),rep(0,90))
  test.genes<-y+c
  Data<-cbind(test.genes,sim.data)
  
#add column and row names
  pre<-"test.gene"
  suf<-seq(1:10)
  out<-paste(pre,suf,sep="")
  prefix<-"gene"
  suffix<-seq(1:1000)
  colnames(Data)<-c(out,paste(prefix,suffix,sep=""))
  rownames(Data)<-paste("sample",seq(1:100),sep="")
  
#apply outlying.samples
  outlying.samples(Data,out)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{outlier}
