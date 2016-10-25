\name{Run_permutation}
\alias{Run_permutation}
\title{Derive importance scores for M permuted data sets.}
\usage{
Run_permutation(X, ntree, mtry,genes.name,M)
}
\arguments{
\item{X}{List object containing expression data for each class, \code{X=list(x_1,x_2, ... )} where \code{x_j} is a \code{(p x n_j)} matrix with rows corresponding to genes and columns to samples. Rows need to be the same across objects, while samples can vary. Missing values are not allowed.}

\item{ntree}{numeric value: number of trees.}

\item{mtry}{numeric value: number of predictors to be sampled at each node.}

\item{genes.name}{vector containing genes name. The order needs to match the rows of \code{x_j}.}

\item{M}{integer: total number of permutations.}

}
\value{
A three dimensional matrix (\code{I} x \code{M} x \code{C}) with \code{I} being the number of total interactions, \code{M} the number of permutations and \code{C} the number of classes. Element \code{(i,j,k)} corresponds to the importance score for interaction \code{i}, permuted data \code{j} and class \code{k}. 
}
\description{
This function computes importance score for \code{M} permuted data sets. Sample labels of target genes are randomly permuted and JRF is implemented. Resulting importance scores can be used to derive an estimate of FDR.  
}
\examples{

 # --- Derive weighted networks via JRF
 
 nclasses=2               # number of data sets / classes
 n1<-n2<-20               # sample size for each data sets
 p<-5                   # number of variables (genes)
 genes.name<-paste("G",seq(1,p),sep="")   # genes name
 perm=1;        # set permutation seed
 
   # --- Generate data sets
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2
 M=5; 
 
   # --- Standardize variables to mean 0 and variance 1
   
  data1 <- t(apply(data1, 1, function(x) { (x - mean(x)) / sd(x) } ))
  data2 <- t(apply(data2, 1, function(x) { (x - mean(x)) / sd(x) } ))
   
   # --- Run JRF and obtain importance score of interactions for each class
   
  out<-Run_permutation(list(data1,data2),mtry=round(sqrt(p-1)),ntree=1000,genes.name,M)

}
\references{
Petralia, F., Song, WM., Tu, Z. and Wang, P., A New Method for Joint Network Analysis Reveals Common and Different Co-Expression Patterns Among Genes and Proteins in Breast Cancer, submitted

A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. \emph{R News} \bold{2}, 18--22.
}
