% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LncFinder.R
\name{lnc_finder}
\alias{lnc_finder}
\title{Long Non-coding RNA Identification}
\usage{
lnc_finder(
  Sequences,
  SS.features = FALSE,
  format = "DNA",
  frequencies.file = "human",
  svm.model = "human",
  parallel.cores = 2
)
}
\arguments{
\item{Sequences}{Unevaluated sequences. Can be a FASTA file loaded by
\code{\link[seqinr]{seqinr-package}} or secondary structure sequences
(Dot-Bracket Notation) obtained from function \code{\link{run_RNAfold}}. If
\code{Sequences} is secondary structure sequences file, parameter \code{format}
should be defined as \code{"SS"}.}

\item{SS.features}{Logical. If \code{SS.features = TRUE}, secondary structure
features will be used.}

\item{format}{String. Define the format of the \code{Sequences}. Can be
\code{"DNA"} or \code{"SS"}. \code{"DNA"} for DNA sequences and \code{"SS"}
for secondary structure sequences.}

\item{frequencies.file}{String or a list obtained from function
\code{\link{make_frequencies}}. Input species name \code{"human"}, \code{"mouse"}
or \code{"wheat"} to use pre-build frequencies files. Or assign a users' own
frequencies file (See function \code{\link{make_frequencies}}).}

\item{svm.model}{String or a svm model obtained from function \code{\link{build_model}}
or \code{\link{svm_tune}}. Input species name \code{"human"}, \code{"mouse"}
or \code{"wheat"} to use pre-build models. Or assign a users' own model (See
function \code{\link{build_model}}).}

\item{parallel.cores}{Integer. The number of cores for parallel computation.
By default the number of cores is \code{2}. Users can set as \code{-1} to run
this function with all cores.}
}
\value{
Returns a data.frame. Including the results of prediction (\code{Pred});
coding potential (\code{Coding.Potential}) and the features. For the details
of the features, please refer to function \code{\link{extract_features}}.
}
\description{
This function is used to predict sequences are non-coding transcripts
or protein-coding transcripts.
}
\details{
Considering that it is time consuming to obtain secondary structure
sequences, users can input nucleotide sequences and predict these sequences
without secondary structure features (Set \code{SS.features} as \code{FALSE}).

Please note that:

\code{SS.features} can improve the performance when the species of unevaluated
sequences is identical to the species of the sequences that used to build the
model.

However, if users are trying to predict sequences with the model trained on
other species, \code{SS.features} may lead to low accuracy.

For the details of \code{frequencies.file}, please refer to function
\code{\link{make_frequencies}}.

For the details of the features, please refer to function
\code{\link{extract_features}}.
}
\section{References}{

Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
LncFinder: an integrated platform for long non-coding RNA identification utilizing
sequence intrinsic composition, structural information, and physicochemical property.
\emph{Briefings in Bioinformatics}, 2019, 20(6):2009-2027.
}

\examples{
\dontrun{
data(demo_DNA.seq)
Seqs <- demo_DNA.seq

### Input one sequence:
OneSeq <- Seqs[1]
result_1 <- lnc_finder(OneSeq, SS.features = FALSE, format = "DNA",
                       frequencies.file = "human", svm.model = "human",
                       parallel.cores = 2)

### Or several sequences:
data(demo_SS.seq)
Seqs <- demo_SS.seq
result_2 <- lnc_finder(Seqs, SS.features = TRUE, format = "SS",
                       frequencies.file = "mouse", svm.model = "mouse",
                       parallel.cores = 2)

### A complete work flow:
### Calculate second structure on Windows OS,
RNAfold.path <- '"E:/Program Files/ViennaRNA/RNAfold.exe"'
SS.seq <- run_RNAfold(Seqs, RNAfold.path = RNAfold.path, parallel.cores = 2)

### Predict the sequences with secondary structure features,
result_2 <- lnc_finder(SS.seq, SS.features = TRUE, format = "SS",
                       frequencies.file = "mouse", svm.model = "mouse",
                       parallel.cores = 2)

### Predict sequences with your own model by assigning a new svm.model and
### frequencies.file to parameters "svm.model" and "frequencies.file"
}
}
\seealso{
\code{\link{build_model}}, \code{\link{make_frequencies}},
         \code{\link{extract_features}}, \code{\link{run_RNAfold}}, \code{\link{read_SS}}.
}
\author{
HAN Siyu
}
