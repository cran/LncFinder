#' A demo of dataset
#'
#' This dataset contains the features of 20 lncRNA seuqences and 20 protein-coding
#' sequences.
#'
#' @format A data frame with 40 rows and 20 variables:
#' \describe{
#'   \item{Label}{the class of the sequences}
#'   \item{ORF.Max.Len}{the length of the longest ORF}
#'   \item{ORF.Max.Cov}{the coverage of the longest ORF}
#'   \item{Seq.lnc.Dist}{Log-Distance.lncRNA}
#'   \item{Seq.pct.Dist}{Log-Distance.protein-coding transcripts}
#'   \item{Seq.Dist.Ratio}{Distance-Ratio.sequence}
#'   \item{Signal.Peak}{Signal as 1/3 position}
#'   \item{SNR}{Signal to noise ratio}
#'   \item{Signal.Min}{the minimum value of the top 10\% power spectrum}
#'   \item{Signal.Q1}{the quantile Q1 of the top 10\% power spectrum}
#'   \item{Signal.Q2}{the quantile Q2 of the top 10\% power spectrum}
#'   \item{Signal.Max}{the maximum value of the top 10\% power spectrum}
#'   \item{Dot_lnc.dist}{Log-Distance.acguD.lncRNA}
#'   \item{Dot_pct.dist}{Log-Distance.acguD.protein-coding transcripts}
#'   \item{Dot_Dist.Ratio}{Distance-Ratio.acguD}
#'   \item{SS.lnc.dist}{Log-Distance.acgu-ACGU.lncRNA}
#'   \item{SS.pct.dist}{Log-Distance.acgu-ACGU.protein-coding transcripts}
#'   \item{SS.Dist.Ratio}{Distance-Ratio.acgu-ACGU}
#'   \item{MFE}{Minimum free energy}
#'   \item{UP.PCT}{Percentage of Unpair-Pair}
#' }
#' @source Sequences are selected from GENCODE.
"demo_dataset"

#' A demo of secondary structure sequences
#'
#' This file contains 10 SS (Secondary Structure) sequences.
#'
#' @format A data frame with 3 rows and 10 variables:
#' \describe{
#' The first row is RNA sequence; the second row is Dot-Bracket Notation of
#' secondary structure sequences; the last row is minimum free energy (MFE).
#' }
#' @source DNA sequences are selected from GENCODE. Secondary struture of each
#' sequence is obtained from progrem "RNAfold".
"demo_SS.seq"

#' A demo of DNA sequences
#'
#' This file contains 10 DNA sequences.
#'
#' @format A list contains 10 DNA sequences.
#' \describe{
#' The sequences are loaded by function \code{\link[seqinr]{read.fasta}}.
#' }
#' @source DNA sequences are selected from GENCODE.
"demo_DNA.seq"
