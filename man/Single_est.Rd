\name{Single_est}
\alias{Single_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Direct estimate in the first step
}
\description{
Utilize the sparsity of the design matrix to estimate partial of betas to lower the dimension of betas.
}
\usage{
Single_est(Counts_Sample, SNP_Pattern, N = nrow(SNP_Pattern), filename = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts_Sample}{
## Original sample data
}
  \item{SNP_Pattern}{
## Original SNP pattern data
}
  \item{N}{
## Row number of SNP pattern
}
  \item{filename}{
## Optional output directory
}
}
\value{
The estimate result from Step I.
}
\author{
Sandeep J. Joseph and Ben Li
}

\examples{
##Load Sample data
data(Sample)
##Load SNP patter data
data(SNP.Pattern)
## Or specify your own data
## Run the first Step
Result2 = Single_est(Sample, SNP.Pattern)
}    
