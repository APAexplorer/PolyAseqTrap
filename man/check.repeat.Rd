% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolyAseqTrap_funclib.R
\name{check.repeat}
\alias{check.repeat}
\title{A check.repeat function}
\usage{
check.repeat(pac.result = NULL, bsgenome = NULL)
}
\arguments{
\item{pac.result}{A data frame of PACs (polyA site clusters), which must include a 'coord' column indicating the genomic coordinate of each representative polyA site.}

\item{bsgenome}{A BSgenome object for reference genome. e.g. 'BSgenome.Hsapiens.UCSC.hg38'.}
}
\value{
A data frame identical to \code{pac.result}, with an additional logical column (e.g., \code{repeat_flag}) indicating whether each representative polyA site is located in an A-rich repeat region.
}
\description{
This function examines whether the representative polyA sites in a PAC data frame
are located within A-rich repeat regions. It uses a reference genome (BSgenome)
to extract the local sequence around each site and assess A-richness and repeat content.
}
