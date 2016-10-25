\name{JRF_network}
\alias{JRF_network}
\title{Compute FDR of importance scores and return class-specific networks.}
\usage{
JRF_network(out.jrf,out.perm,TH)
}
\arguments{
\item{out.jrf}{output object from function JRF.}

\item{out.perm}{output object from function Run_permutation.}

\item{TH}{Threshold for FDR.}
}

\value{
\code{out} list object containing the estimated gene-gene interactions for each class. 
}
\description{
This function computes FDR of importance scores and returns class-specific networks.  
}
\examples{

 # --- Derive weighted networks via JRF
 
 nclasses=2             # number of data sets / classes
 n1<-n2<-20             # sample size for each data sets
 p<-5                   # number of variables (genes)
 genes.name<-paste("G",seq(1,p),sep="")   # genes name
 M=5;                   # total number of permutations
 fdr=.001;              # fdr threshold
 
   # --- Generate data sets
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2
 data1[1,]<-2*data1[2,]     # variable 1 and 2 interact under class 1
  
  # --- Standardize variables to mean 0 and variance 1
   
  data1 <- t(apply(data1, 1, function(x) { (x - mean(x)) / sd(x) } ))
  data2 <- t(apply(data2, 1, function(x) { (x - mean(x)) / sd(x) } ))
   
   # --- Run JRF and obtain importance score of interactions for each class
  
  out<-JRF(list(data1,data2),mtry=round(sqrt(p-1)),ntree=1000,genes.name)
  
  out.perm<-Run_permutation(list(data1,data2),mtry=round(sqrt(p-1)),ntree=1000,genes.name,M)

  final.net<-JRF_network(out,out.perm,fdr)

}
\references{
Petralia, F., Song, WM., Tu, Z. and Wang, P., A New Method for Joint Network Analysis Reveals Common and Different Co-Expression Patterns Among Genes and Proteins in Breast Cancer, submitted

A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. \emph{R News} \bold{2}, 18--22.

Xie, Y., Pan, W. and Khodursky, A.B., 2005. A note on using permutation-based false discovery rate estimates to compare different analysis methods for microarray data. \emph{Bioinformatics}, \bold{21}(23), pp.4280-4288.
}
