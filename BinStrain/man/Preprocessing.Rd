\name{Preprocessing}
\alias{Preprocessing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preprocess the original data.
}
\description{
Convert the original data information into a design matrix. Deal with missing data.
}
\usage{
Preprocessing(Counts_Sample, SNP_Pattern, N = nrow(SNP_Pattern))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts_Sample}{
## Original Sample data
}
  \item{SNP_Pattern}{
## Original SNP pattern data
}
  \item{N}{
## row number of SNP pattern
}
}
\value{
Return a list containing the information from original data with a better format
}
\author{
Sandeep J. Joseph and Ben Li
}

\examples{
## First loading our example data or specify your own data
## Run the function
Sample2 = Preprocessing(Sample, SNP.Pattern)
}

