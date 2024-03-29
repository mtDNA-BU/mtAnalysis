% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.mtDNAaaf.R
\name{plot.mtDNAaaf}
\alias{plot.mtDNAaaf}
\title{plot.mtDNAaaf function}
\usage{
\method{plot}{mtDNAaaf}(x, col = "blue", pch = ".", cex = 0.2, xlab = "", ylab = "", ...)
}
\arguments{
\item{x}{a numeric matrix (16569 x N).
Rows correspond to loci and columns correspond to subjects.
It contains subject ID as the column names,
and the AAFs of all 16569 mtDNA loci for each subject.
It is generated from mtAAF function.}

\item{col}{color to be used in the plot.}

\item{pch}{plotting "character", i.e., symbol to use.}

\item{cex}{expansion factor for symbols used in the plot.}

\item{xlab}{a label for the x axis.}

\item{ylab}{a label for the y axis.}

\item{...}{arguments to be passed to plot function.}
}
\description{
Produce a Diagnostic Scatter Plot to Visualize AAFs Generated by
mtAAF Function
}
\examples{

\dontrun{
## Read input data
allele_file <- "allele.csv"
freq_file   <- "freq.csv"

allele <- as.matrix( read.csv(file=allele_file, sep=",", header=FALSE))
freq   <- as.matrix( read.csv(file=freq_file, sep=",", header=FALSE))

aaf=mtAAF(allele, freq)
plot(aaf)
}

}
