% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LncFinder.R
\name{compute_kmer}
\alias{compute_kmer}
\title{Compute \emph{k}-mer Features}
\usage{
compute_kmer(
  Sequences,
  label = NULL,
  k = 1:5,
  step = 1,
  freq = TRUE,
  improved.mode = FALSE,
  alphabet = c("a", "c", "g", "t"),
  on.ORF = FALSE,
  auto.full = FALSE,
  parallel.cores = 2
)
}
\arguments{
\item{Sequences}{A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
\code{\link[seqinr]{seqinr-package}}.}

\item{label}{Optional. String. Indicate the label of the sequences such as
"NonCoding", "Coding".}

\item{k}{An integer that indicates the sliding window size. (Default: \code{1:5})}

\item{step}{Integer defaulting to \code{1} for the window step.}

\item{freq}{Logical. If TRUE, the frequencies of different patterns are returned
instead of counts. (Default: \code{TRUE})}

\item{improved.mode}{Logical. If TRUE, the frequencies will be normalized using
the method proposed by PLEK (Li et al. 2014).
Ignored if \code{freq = FALSE}. (Default: \code{FALSE})}

\item{alphabet}{A vector of single characters that specify the different character
of the sequence. (Default: \code{alphabet = c("a", "c", "g", "t")})}

\item{on.ORF}{Logical. If \code{TRUE}, the \emph{k}-mer frequencies will be calculated on
the longest ORF region. NOTE: If \code{TRUE}, the sequences have to be DNA.
(Default: \code{FALSE})}

\item{auto.full}{Logical. When \code{on.ORF = TRUE} but no ORF can be found,
if \code{auto.full = TRUE}, the \emph{k}-mer
frequencies will be calculated on the full sequence automatically;
if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded.
Ignored when \code{on.ORF = FALSE}. (Default: \code{FALSE})}

\item{parallel.cores}{Integer. The number of cores for parallel computation.
By default the number of cores is \code{2}. Users can set as \code{-1} to run
this function with all cores.}
}
\value{
A dataframe.
}
\description{
This function can calculate the \emph{k}-mer frequencies of the sequences.
}
\details{
This function can extract \emph{k}-mer features. \code{k} and \code{step} can be customized.
The count (\code{freq = FALSE}) or frequencies (\code{freq = TRUE}) of different patterns can be returned.
If \code{freq = TRUE}, \code{improved.mode} is available. The improved mode is proposed by method PLEK.
(Ref: Li et al. 2014)
}
\examples{
\dontrun{
data(demo_DNA.seq)
Seqs <- demo_DNA.seq

kmer_res1 <- compute_kmer(Seqs, k = 1:5, step = 1, freq = TRUE, improved.mode = FALSE)

kmer_res2 <- compute_kmer(Seqs, k = 1:5, step = 3, freq = TRUE,
                          improved.mode = TRUE, on.ORF = TRUE, auto.full = TRUE)
}

}
\author{
HAN Siyu
}
