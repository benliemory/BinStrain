\name{Whole_est}
\alias{Whole_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Second Step of Estimate
}
\description{
Using quadratic optimization method to estimate all the betas.
}
\usage{
Whole_est(results, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
The results from the estimate in Step I.
}
  \item{data}{
The data after preprocessing.
}
}
\value{
Return the estimate after step II.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Sandeep J. Joseph and Ben Li
}
\examples{
##Load the example data or specify your own data
data(Sample)
data(SNP.Pattern)
##Preprocess the sample data
Sample2 = Preprocessing(Sample, SNP.Pattern)
##Step I estimate
Result2 = Single_est(Sample, SNP.Pattern)
##Step II estimate
betas = Whole_est(Result2,Sample2)
}
