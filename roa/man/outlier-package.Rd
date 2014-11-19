\name{roa}
\alias{roa}
\docType{package}
\title{Robust Outlier Analysis}
\description{
This package has the capabilities to determine outliers from a variety of biological data features, such as microarray gene expression, RNA-Seq, methylation, and copy number. Outliers are identified using previously developed univariate outlier detection methods: the Cancer Profile Analysis (COPA), Outlier Sum Statistic (OSS), Gene Tissue Index (GTI), and a simple variance approach. Outliers are also determined by a novel multivariate approach. Both the univariate and multivariate approaches are used in concert with either a change point model or a user-defined cut-off to determine the threshold above, or below, which outliers are identified. Univariate outliers can be intersected to find outliers in common across methods.
}
\details{
Use uv.outlier to find univariate outliers in one type of sample, for example all disease samples. Use uv.outlier2 to determine univariate outliers in disease samples versus normal samples. Use mv.outlier to identify multivariate outliers. To identify outliers using only one univariate approach use, copa, gti, oss, or variance.
}
\author{
Kylie Ainslie and Jeanne Kowalski
}
\references{
Mpindi, J.P. et al. 2011. GTI: A novel algorithm for identifying outlier gene expression profiles from integrated microarray datasets.PLoS One 6(2): e17259.
Ross, G.J. Parametric and nonparametric sequential change detection in R: the cpm package. Journal of Statistical Software  (forthcoming).
Tibshirani, R. and T. Hastie. 2007. Outlier sums for differential gene expression analysis. Biostatistics 8(1): 2-8.
Tomlins, S.A. et al. 2005. Recurrent fusion of TMPRSS2 and ETS transcription factor genes in prostate cancer. Science 310: 644-648

}

\keyword{ package }
