\name{BinStrain-package}
\alias{BinStrain-package}
\alias{BinStrain}
\docType{package}
\title{
BinStrain Overview
}
\description{
Description: Uses a binomial model to describe the observed alternative allele (SNP) in the targeted genomic region (or the whole genome). A well-established two step method is used to estimate the proportion of serovar strain-specific SNPs present in pure or mixed infection of genotypes. First step is a direct estimate by using the sparsity of the design matrix. Quadratic optimization method is involved in the second step. Package quadprog is needed.
}
\details{
\tabular{ll}{
Package: \tab BinStrain\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-28\cr
License: \tab GPL-2\cr
}
After loading your data with appropriate format, using the function "Main.function" can directly generate a text file named "betas.txt" which contains all the estimate for the parameter beta in our model.The first column of the "betas.txt" file represents the assigned consecutive numbers used to represent the reference genomes used for creating the SNP pattern file. If there are 2 or more highly similar strains with out any differentiating SNPs between them, then beta.txt file will have beta estimate for only one strain and for the rest of the similar strains will have 'NA' as estimate.
}
\author{
Sandeep J. Joseph and Ben Li

Maintainer: Sandeep J. Joseph <sandeepjoseph@emory.edu> and Ben Li <ben.li@emory.edu>
}

\keyword{ package }
