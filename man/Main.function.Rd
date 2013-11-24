\name{Main.function}
\alias{Main.function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to estimate the parameters in our model
}
\description{
Sample data and SNP pattern data are needed. After loading the data with an appropriate format, this function can generate the results directly.
}
\usage{
Main.function(Sample, SNP.Pattern)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sample}{
The data for sample. The data should follow the format of our sample data. See data(Sample) more detailed format.    
}
  \item{SNP.Pattern}{
The data for SNP pattern. The data should follow the format of our SNP pattern data.See data(SNP.Pattern) for more detailed format. 
}
}

\value{
No values returned. A text file for all the betas value will be generated in the current folder.
}
\author{
Sandeep J. Joseph and Ben Li
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##Load sample data
data(Sample)
##Load SNP pattern data
data(SNP.Pattern)
## Or specify your own data
##Run the main function
Main.function(Sample,SNP.Pattern) 
}

