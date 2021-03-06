% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LncFinder.R
\name{compute_LogDistance}
\alias{compute_LogDistance}
\title{Compute Logarithm Distance}
\usage{
compute_LogDistance(
  Sequences,
  label = NULL,
  referFreq,
  k = 6,
  step = 1,
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

\item{referFreq}{a list obtained from function \code{\link{make_referFreq}}.}

\item{k}{An integer that indicates the sliding window size. (Default: \code{6})}

\item{step}{Integer defaulting to \code{1} for the window step.}

\item{alphabet}{A vector of single characters that specify the different character
of the sequence. (Default: \code{alphabet = c("a", "c", "g", "t")})}

\item{on.ORF}{Logical. If \code{TRUE}, Logarithm Distance will be calculated on
the longest ORF region. NOTE: If \code{TRUE}, the input has to be DNA sequences.
(Default: \code{FALSE})}

\item{auto.full}{Logical. When \code{on.ORF = TRUE} but no ORF can be found,
if \code{auto.full = TRUE}, Logarithm Distance will be calculated on full sequences automatically;
if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded. Ignored when \code{on.ORF = FALSE}.
(Default: \code{FALSE})}

\item{parallel.cores}{Integer. The number of cores for parallel computation.
By default the number of cores is \code{2}. Users can set as \code{-1} to run
this function with all cores.}
}
\value{
A dataframe.
}
\description{
This function can compute Logarithm Distance proposed by method LncFinder
(Han et al. 2018). Logarithm Distance can be calculated on full sequence or the longest ORF
region. The step and \emph{k} of the sliding window can also be customized.
}
\details{
This function can compute Logarithm Distance proposed by LncFinder (HAN et al. 2018).
In LncFinder, two schemes are provided to calculate Logarithm Distance:
1) \code{step = 3} and \code{k = 6} on the longest ORF region;
2) \code{step = 1} and \code{k = 6} on full sequence.
Method LncFinder uses scheme 1 to extract Logarithm Distance features.
Using this function \code{compute_EucDistance}, both \code{step}, \code{k},
and calculated region (full sequence or ORF)
can be customized to maximize its availability.
}
\section{References}{

Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
LncFinder: an integrated platform for long non-coding RNA identification utilizing
sequence intrinsic composition, structural information, and physicochemical property.
\emph{Briefings in Bioinformatics}, 2019, 20(6):2009-2027.
}

\examples{
\dontrun{
Seqs <- seqinr::read.fasta(file =
"http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")

referFreq <- make_referFreq(cds.seq = Seqs, lncRNA.seq = Seqs, k = 6, step = 3,
                            alphabet = c("a", "c", "g", "t"), on.orf = TRUE,
                            ignore.illegal = TRUE)

data(demo_DNA.seq)
Sequences <- demo_DNA.seq

LogDistance <- compute_LogDistance(Sequences, label = "NonCoding", referFreq = referFreq,
                                   k = 6, step = 3, alphabet = c("a", "c", "g", "t"),
                                   on.ORF = TRUE, auto.full = TRUE, parallel.cores = 2)
}

}
\seealso{
\code{\link{make_referFreq}},
         \code{\link{compute_EucDistance}},
         \code{\link{compute_hexamerScore}}.
}
\author{
HAN Siyu
}
