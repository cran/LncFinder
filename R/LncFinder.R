################  Internal Functions  ################

get.Seq_UP <- function(OneSeq) {
        UP.seq <- seqinr::s2c(OneSeq[[2]])
        UP.seq[UP.seq %in% "(" | UP.seq %in% ")"] <- "P"
        UP.seq[UP.seq %in% "."] <- "U"
        UP.seq
}

get.Seq_acguACGU <- function(OneSeq) {
        ACGU.seq   <- seqinr::s2c(OneSeq[[1]])
        Stem.index <- seqinr::s2c(OneSeq[[2]]) %in% "(" | seqinr::s2c(OneSeq[[2]]) %in% ")"
        ACGU.seq[Stem.index & ACGU.seq %in% "a"] <- "A"
        ACGU.seq[Stem.index & ACGU.seq %in% "c"] <- "C"
        ACGU.seq[Stem.index & ACGU.seq %in% "g"] <- "G"
        ACGU.seq[Stem.index & ACGU.seq %in% "u"] <- "U"
        ACGU.seq
}

get.Seq_acguD <- function(OneSeq) {
        acguD.seq <- seqinr::s2c(OneSeq[[1]])
        Dot.index <- seqinr::s2c(OneSeq[[2]]) %in% "."
        acguD.seq[Dot.index] <- "D"
        acguD.seq
}

get.Seq_DNA <- function(OneSeq) {
        DNA.seq <- seqinr::s2c(OneSeq[[1]])
        DNA.seq[DNA.seq %in% "u"] <- "t"
        DNA.seq
}

find_ORF <- function(OneSeq) {
        OneSeq <- unlist(seqinr::getSequence(OneSeq, as.string = TRUE))
        start_pos <- unlist(gregexpr("atg", OneSeq, ignore.case = TRUE))
        if(sum(start_pos) == -1) {
                max_len <- 0
                max_cov <- 0
                orf_seq <- NA
        } else {
                stop_pos <- sapply(c("taa", "tag", "tga"), function(x){
                        pos <- unlist(gregexpr(x, OneSeq, ignore.case = TRUE))
                })
                stop_pos <- sort(unlist(stop_pos, use.names = FALSE))
                seq_length <- nchar(OneSeq)
                if(all(stop_pos == -1)){
                        orf_start <- min(start_pos)
                        orf_stop  <- seq_length - ((seq_length - orf_start + 1) %% 3) - 2
                        max_len   <- orf_stop - orf_start + 3
                        max_cov   <- max_len / seq_length
                        orf_seq   <- substr(OneSeq, start = orf_start, stop = orf_stop)
                } else {
                        orf_starts  <- c()
                        orf_stops   <- c()
                        orf_lengths <- c()
                        for(i in start_pos){
                                orf_flag <- 0
                                for(j in stop_pos){
                                        if(j < i) next
                                        diff_mod    <- (j - i) %% 3
                                        if(diff_mod != 0) next
                                        orf_starts  <- c(orf_starts, i)
                                        orf_stops   <- c(orf_stops, j)
                                        orf_lengths <- c(orf_lengths, j - i + 3)
                                        orf_flag    <- 1
                                        break
                                }
                                if(orf_flag == 0){
                                        orf_starts   <- c(orf_starts, i)
                                        orf_stop_tmp <- seq_length - ((seq_length - i + 1) %% 3) - 2
                                        orf_stops    <- c(orf_stops, orf_stop_tmp)
                                        orf_lengths  <- c(orf_lengths, orf_stop_tmp - i + 3)
                                }
                        }
                        max_len   <- max(orf_lengths)
                        max_cov   <- max_len / seq_length
                        max_index <- which(orf_lengths == max_len)
                        orf_start <- orf_starts[max_index]
                        orf_stop  <- orf_stops[max_index]
                        orf_seq   <- substr(OneSeq, start = orf_start, stop = orf_stop + 2)
                }

        }
        ORF_info <- list(ORF.Max.Len = max_len, ORF.Max.Cov = max_cov, ORF.Max.Seq = orf_seq)
        ORF_info
}

get.freq <- function(seqs, alphabet, wordsize, step, freq, ignore.illegal = TRUE) {

        if(ignore.illegal) {
                for(i in 1:length(seqs)) {
                        if(!all(seqs[[i]] == "a" | seqs[[i]] == "c" | seqs[[i]] == "g" | seqs[[i]] == "t")) next
                        freq_tmp  <- seqinr::count(seqs[[i]], wordsize = wordsize,
                                                   by = step, freq = freq, alphabet = alphabet)
                        if(i == 1) freq.res <- freq_tmp else freq.res <- freq.res + freq_tmp
                        if(i %% 5000 == 0) message("  ", i, " sequences completed.")
                }
        } else {
                for(i in 1:length(seqs)) {
                        freq_tmp  <- seqinr::count(seqs[[i]], wordsize = wordsize,
                                                   by = step, freq = freq, alphabet = alphabet)
                        if(i == 1) freq.res <- freq_tmp else freq.res <- freq.res + freq_tmp
                        if(i %% 5000 == 0) message("  ", i, " sequences completed.")
                }
        }

        freq.out <- freq.res / sum(freq.res)
        freq.out
}

LogDist.DNA  <- function(OneSeq, hexamer.lnc, hexamer.cds, wordsize, step) {
        orf.info  <- find_ORF(OneSeq)

        count6  <- seqinr::count(seqinr::s2c(orf.info[[3]]), wordsize = wordsize,
                                 by = step, freq = FALSE)

        if(sum(count6) < 3) {
                count6  <- seqinr::count(OneSeq,  wordsize = wordsize,
                                         by = step, freq = FALSE)
        }

        freq6     <- count6 / sum(count6)
        lnc.ratio <- freq6[hexamer.lnc != 0] / hexamer.lnc[hexamer.lnc != 0]
        pct.ratio <- freq6[hexamer.cds != 0] / hexamer.cds[hexamer.cds != 0]
        lnc.dist  <- (sum(log(lnc.ratio[lnc.ratio != 0])) - length(which(lnc.ratio == 0)) + length(which(hexamer.lnc == 0))) / sum(count6)
        pct.dist  <- (sum(log(pct.ratio[pct.ratio != 0])) - length(which(pct.ratio == 0)) + length(which(hexamer.cds == 0))) / sum(count6)
        Ratio     <- lnc.dist / pct.dist
        distance  <- c(ORF.Max.Len = orf.info[[1]], ORF.Max.Cov = orf.info[[2]],
                       Seq.lnc.Dist = lnc.dist, Seq.pct.Dist = pct.dist, Seq.Dist.Ratio = Ratio)
}

LogDist.acguD  <- function(OneSeq, hexamer.lnc, hexamer.cds, wordsize, step) {
        count6    <- seqinr::count(OneSeq, wordsize = wordsize, freq = FALSE, by = step,
                                   alphabet = c("D", "a", "c", "g", "u"))
        freq6     <- count6 / sum(count6)
        lnc.ratio <- freq6[hexamer.lnc != 0] / hexamer.lnc[hexamer.lnc != 0]
        pct.ratio <- freq6[hexamer.cds != 0] / hexamer.cds[hexamer.cds != 0]
        lnc.dist  <- (sum(log(lnc.ratio[lnc.ratio != 0])) - length(which(lnc.ratio == 0)) + length(which(hexamer.lnc == 0))) / sum(count6)
        pct.dist  <- (sum(log(pct.ratio[pct.ratio != 0])) - length(which(pct.ratio == 0)) + length(which(hexamer.cds == 0))) / sum(count6)
        Ratio     <- lnc.dist / pct.dist
        distance  <- c(Dot_lnc.dist = lnc.dist, Dot_pct.dist = pct.dist, Dot_Dist.Ratio = Ratio)
}

LogDist.acguACGU  <- function(OneSeq, hexamer.lnc, hexamer.cds, wordsize, step) {
        count6    <- seqinr::count(OneSeq, wordsize = wordsize, freq = FALSE, by = step,
                                   alphabet = c("A", "C", "G", "U", "a", "c", "g", "u"))
        freq6     <- count6 / sum(count6)
        lnc.ratio <- freq6[hexamer.lnc != 0] / hexamer.lnc[hexamer.lnc != 0]
        pct.ratio <- freq6[hexamer.cds != 0] / hexamer.cds[hexamer.cds != 0]
        lnc.dist  <- (sum(log(lnc.ratio[lnc.ratio != 0])) - length(which(lnc.ratio == 0)) + length(which(hexamer.lnc == 0))) / sum(count6)
        pct.dist  <- (sum(log(pct.ratio[pct.ratio != 0])) - length(which(pct.ratio == 0)) + length(which(hexamer.cds == 0))) / sum(count6)
        Ratio     <- lnc.dist / pct.dist
        distance  <- c(SS.lnc.dist = lnc.dist, SS.pct.dist = pct.dist, SS.Dist.Ratio = Ratio)
}

EIIP.DFT <- function(OneSeq) {

        OneSeq <- seqinr::getSequence(OneSeq)
        OneSeq[!(OneSeq %in% "a" | OneSeq %in% "c" | OneSeq %in% "g" | OneSeq %in% "t")] <- 0
        OneSeq[OneSeq %in% "a"] <- 0.1260
        OneSeq[OneSeq %in% "c"] <- 0.1340
        OneSeq[OneSeq %in% "g"] <- 0.0806
        OneSeq[OneSeq %in% "t"] <- 0.1335

        numeric.seq  <- as.numeric(OneSeq)
        fourier.seq  <- stats::fft(numeric.seq)
        seq.spectrum <- (abs(fourier.seq)) ^ 2
        power.length <- length(seq.spectrum)
        k <- round(power.length / 3)
        DFT.max <- max(seq.spectrum[k-2], seq.spectrum[k-1], seq.spectrum[k],
                       seq.spectrum[k+1], seq.spectrum[k+2])
        total.power  <- sum(seq.spectrum) / power.length
        SNR.res <- c(Signal.Peak = DFT.max, SNR = (DFT.max / total.power))

        seq.spectrum <- seq.spectrum[-c(1, power.length)]
        seq.spectrum <- unique(seq.spectrum)
        ordered.pow  <- sort(seq.spectrum, decreasing = TRUE)
        pow.percent  <- ordered.pow[1:round(length(ordered.pow) * 0.1)]
        Quantile.res <- stats::quantile(pow.percent)[-4]
        names(Quantile.res) <- c("Signal.Min", "Signal.Q1", "Signal.Q2", "Signal.Max")
        EIIP.res <- c(SNR.res, Quantile.res)

        return(EIIP.res)
}

secondary_seq <- function(OneSeq, info, RNAfold.path, index){
        seq.string <- unlist(seqinr::getSequence(OneSeq, TRUE))

        cat(index, "of", info, "length:", nchar(seq.string), "nt")

        seq.ss <- system(RNAfold.path, intern = TRUE, input = seq.string)
        index <<- index + 1
        seq.ss[3] <- as.numeric(substr(seq.ss[2], nchar(seq.string) + 3, nchar(seq.ss[2]) - 1))
        seq.ss[2] <- substr(seq.ss[2], 1, nchar(seq.string))
        seq.ss
}

################  Main Functions

###### max_orf ######

#' Find the longest ORF
#' @description This function can find the longest ORF in one sequence.
#'
#' @param OneSeq Is one sequence. Can be a FASTA file read by package "seqinr"
#' (\code{\link[seqinr]{seqinr-package}}) or just a string.
#'
#' @param reverse.strand Logical. Whether find ORF on the reverse strand.
#'
#' @return Returns a list which consists ORF region (\code{ORF.Max.Seq}), length
#' (\code{ORF.Max.Len}) and coverage (\code{ORF.Max.Cov}) of the longest ORF.
#' @author Han Siyu
#' @details This function can extract the longest ORF of one sequence. It returns
#' the region, length and coverage of the longest ORF. Coverage is the the ratio
#' of ORF to transcript length. If \code{reverse.strand = TRUE}, ORF will be
#' found on both the forward and reverse strand.
#' @importFrom seqinr getSequence
#' @importFrom seqinr comp
#' @examples
#' ### For one sequence:
#' OneSeq <- c("cccatgcccagctagtaagcttagcc")
#' max_orf_1 <- max_orf(OneSeq, reverse.strand = TRUE)
#'
#' ### For a FASTA file contains several sequences:
#' \dontrun{
#' ### Use "read.fasta" function of package "seqinr" to read a FASTA file:
#' Seqs <- read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#' }
#'
#' ### Or just try to use our data "demo_DNA.seq"
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' ### Use apply function to find the longest ORF:
#' max_orf_2 <- sapply(Seqs, max_orf, reverse.strand = FALSE)
#' @export

max_orf <- function(OneSeq, reverse.strand = FALSE) {

        orf.info <- find_ORF(OneSeq)
        if(reverse.strand) {
                OneSeq <- unlist(seqinr::getSequence(OneSeq))
                Seq.reverse <- rev(seqinr::comp(OneSeq))
                orf.reverse <- find_ORF(Seq.reverse)
                orf.info    <- list(ORF.Forward = orf.info,
                                    ORF.Reverse = orf.reverse)
        }
        orf.info
}

###### run_RNAfold ######

#' Obtain the Secondary Structure Sequences Using RNAfold
#' @description This function can compute secondary structure sequences. The tool
#' "RNAfold" of software "ViennaRNA" is required for this function.
#'
#' @param Sequences A FASTA file loaded by function \code{\link{read.fasta}} of
#' package "seqinr" (\code{\link[seqinr]{seqinr-package}}).
#'
#' @param RNAfold.path String. Indicate the path of the program "RNAfold". By
#' default is \code{"RNAfold"} for UNIX/Linux system. (See details.)
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}, users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return Returns data.frame. The first row is RNA sequence; the second row is
#' Dot-Bracket Notation of secondary structure sequences; the last row is minimum
#' free energy (MFE).
#' @author Han Siyu
#' @details This function is used to compute secondary structure. The output of
#' this function can be used in function \code{\link{make_frequencies}},
#' \code{\link{extract_features}}, \code{\link{build_model}} and
#' \code{\link{lnc_finder}} when parameter \code{SS.features} is set as \code{TRUE}.
#'
#' This function depends on the program "RNAfold" of software "ViennaRNA".
#' (\url{http://www.tbi.univie.ac.at/RNA/index.html})
#'
#' Parameter \code{RNAfold.path} can be simply defined as \code{"RNAfold"} as
#' default when the OS is UNIX/Linux. However, for some OS, such as Windows, users
#' need to specify the \code{RNAfold.path} if the path of "RNAfold" haven't been
#' added in environment variables.
#'
#' This function can print the related information when the OS is UNIX/Linux,
#' such as:
#'
#' \code{"25 of 100, length: 695 nt"},
#'
#' which means around 100 sequences are assigned to this node and the program is
#' computing the 25th sequence. The length of this sequence is 695nt.
#'
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @examples
#' \dontrun{
#' ### For a FASTA file contains several sequences,
#' ### Use "read.fasta" function of package "seqinr" to read a FASTA file:
#' Seqs <- read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#'
#' ### Or just try to use our data "demo_DNA.seq"
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' ### Windows:
#' RNAfold.path <- '"E:/Program Files/ViennaRNA/RNAfold.exe"'
#' SS.seq_1 <- run_RNAfold(Seqs, RNAfold.path = RNAfold.path, parallel.cores = 2)
#'
#' ### For UNIX/Linux, "RNAfold.path" can be just defined as "RNAfold" as default:
#' SS.seq_2 <- run_RNAfold(Seqs, RNAfold.path = "RNAfold", parallel.cores = 2)
#' }
#' @export

run_RNAfold <- function(Sequences, RNAfold.path = "RNAfold",
                        parallel.cores = 2){

        Sequences <- Sequences[which(lengths(Sequences) < 30000)]

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!")

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores, outfile = "")

        message("Sequences Number: ", length(Sequences))
        message("Processing.")
        index <- 1
        info <- paste(round(length(Sequences) / parallel.cores), ",", sep = "")
        parallel::clusterExport(cl, varlist = c("info", "index"), envir = environment())
        sec.seq <- parallel::parSapply(cl, Sequences, secondary_seq, info = info,
                                       RNAfold.path = RNAfold.path, index = index)
        parallel::stopCluster(cl)
        sec.seq <- data.frame(sec.seq, stringsAsFactors = FALSE)
        message("Completed.")
        sec.seq
}

###### make_frequencies ######

#' Make the frequencies file with your own dataset
#' @description This function is used to calculate the frequencies of CDs and
#' secondary structure sequences. This function is useful when users are trying
#' to build their own model.
#'
#' @param cds.seq Coding sequences (mRNA without UTRs). Can be a FASTA file loaded
#' by package "seqinr" (\code{\link[seqinr]{seqinr-package}}) or secondary structure
#' sequences (Dot-Bracket Notation) obtained form function \code{\link{run_RNAfold}}.
#' CDs are used to calculate hexamer frequencies of nucleotide sequences,thus
#' secondary structure is not needed. Parameter \code{cds.format} should be
#' \code{"SS"} when input is secondary structure sequences. (See details for
#' more information.)
#'
#' @param mRNA.seq mRNA sequences with Dot-Bracket Notation. The secondary
#' structure sequences can be obtained from function \code{\link{run_RNAfold}}.
#' mRNA sequences are used to calculate the frequencies of acgu-ACGU and a acguD
#' (see details), thus, mRNA sequences are required only when \code{SS.features = TRUE}.
#'
#' @param lncRNA.seq Long non-coding RNA sequences. Can be a FASTA file loaded by
#' package "seqinr" (\code{\link[seqinr]{seqinr-package}}) or secondary structure
#' sequences (Dot-Bracket Notation) obtained from function \code{\link{run_RNAfold}}.
#' If \code{SS.features = TRUE}, \code{lncRNA.seq} must be RNA sequences with
#' secondary structure sequences and parameter \code{lnc.format} should be defined
#' as \code{"SS"}.
#'
#' @param SS.features Logical. If \code{SS.features = TRUE}, frequencies of secondary
#' structure will also be calculated and the model can be built with secondary
#' structure features. In this case, \code{mRNA.seq} and \code{lncRNA.seq} should
#' be secondary structure sequences.
#'
#' @param cds.format String. Define the format of the sequences of \code{cds.seq}.
#' Can be \code{"DNA"} or \code{"SS"}. \code{"DNA"} for DNA sequences and \code{"SS"}
#' for secondary structure sequences.
#'
#' @param lnc.format String. Define the format of lncRNAs (\code{lncRNA.seq}).
#' Can be \code{"DNA"} or \code{"SS"}. \code{"DNA"} for DNA sequences and \code{"SS"}
#' for secondary structure sequences. This parameter must be defined as \code{"SS"}
#' when \code{SS.features = TURE}.
#'
#' @param check.cds Logical. Incomplete CDs can lead to a false shift and a
#' inaccurate hexamer frequencies. When \code{check.cds = TRUE}, hexamer frequencies
#' will be calculated on the longest ORF. This parameter is strongly recommended to
#' set as \code{TRUE} when mRNA is used as CDs.
#'
#' @param ignore.illegal Logical. If \code{TRUE}, the sequences with non-nucleotide
#' characters (nucleotide characters: "a", "c", "g", "t") will be ignored when
#' calculating hexamer frequencies.
#'
#' @return Returns a list which consists the frequencies of protein-coding sequences
#' and non-coding sequences.
#' @author Han Siyu
#'
#' @details This funtion is used to make frequencies file. This file is needed
#' when users are trying to build their own model.
#'
#' In order to achieve high accuracy, mRNA should not be regarded as CDs and assigned
#' to parameter \code{cds.seq}. However, CDs of some species may be insufficient
#' for calculating frequencies, and mRNAs can be regared as CDs with parameter
#' \code{check.cds = TRUE}. In this case, hexamer frequencies will be calculated
#' on ORF region.
#'
#' Considering that it is time consuming to obtain secondary structure sequences,
#' users can only provide nucleotide sequences and build a model without secondary
#' structure features (\code{SS.features = FALSE}). If users want to build a model
#' with secondary structure features, parameter \code{SS.features} should be set
#' as \code{TRUE}. At the same time, the format of the sequences of \code{mRNA.seq}
#' and \code{lnc.seq} should be secondary structure sequences (Dot-Bracket Notation).
#' Secondary structure sequences can be obtained by function \code{\link{run_RNAfold}}.
#'
#' Please note that:
#'
#' SS.features can improve the performance when the species of unevaluated sequences
#' is identical to the species of the sequences that used to build the model.
#'
#' However, if users are trying to predict sequences with the model trained on
#' other species, SS.features may lead to low accuracy.
#'
#' The frequencies file consists three groups: Hexamer Frequencies; acgu-ACGU
#' Frequencies and acguD Frequencies.
#'
#' Hexamer Frequencies are calculated on the original nucleotide sequences by
#' employing \emph{k}-mer scheme (\emph{k} = 6), and the sliding window will slide
#' 3nt each step.
#'
#' For any secondary structure sequences (Dot-Bracket Notation), if one position
#' is a dot, the corresponding nucleotide of the RNA sequence will be replaced
#' with character "D". acguD Frequencies are the \emph{k}-mer frequencies
#' (\emph{k} = 4) calculated on this new sequences.
#'
#' Similarly, for any secondary structure sequences (Dot-Bracket Notation), if
#' one position is "(" or ")", the corresponding nucleotide of the RNA sequence
#' will be replaced with upper case ("A", "C", "G", "U").
#'
#' A brief example,
#'
#' DNA Sequence:\code{          5'-   t  a  c  a  g  t  t  a  t  g   -3'}
#'
#' RNA Sequence:\code{          5'-   u  a  c  a  g  u  u  a  u  g   -3'}
#'
#' Dot-Bracket Sequence:\code{     5'-   .  .  .  .  (  (  (  (  (  (   -3'}
#'
#' acguD Sequence:\code{         \{     D, D, D, D, g, u, u, a, u, g   \}}
#'
#' acgu-ACGU Sequence:\code{    \{     u, a, c, a, G, U, U, A, U, G   \}}
#'
#' @section References:
#' HAN Siyu, LIANG Yanchun and LI Ying*.
#' LncFinder: Long Non-coding RNA Identification Tool Based on Features of
#' Sequence Intrinsic Composion, Secondary Structure and EIIP Values. (2017)
#' (\emph{Submitted})
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#' @seealso \code{\link{run_RNAfold}}, \code{\link{build_model}},
#'          \code{\link{extract_features}}.
#' @examples
#' ### Only some brief examples:
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' \dontrun{
#' ### Obtain the secondary structure sequences (Windows OS):
#' RNAfold.path <- '"E:/Program Files/ViennaRNA/RNAfold.exe"'
#' SS.seq <- run_RNAfold(Seqs, RNAfold.path = RNAfold.path, parallel.cores = 2)
#'
#' ### Make frequencies file with secondary strucutre features,
#' my_file_1 <- make_frequencies(cds.seq = SS.seq, mRNA.seq = SS.seq,
#'                             lncRNA.seq = SS.seq, SS.features = TRUE,
#'                             cds.format = "SS", lnc.format = "SS",
#'                             check.cds = TRUE, ignore.illegal = FALSE)
#' }
#'
#' ### Make frequencies file without secondary strucutre features,
#' my_file_2 <- make_frequencies(cds.seq = Seqs, lncRNA.seq = Seqs,
#'                               SS.features = FALSE, cds.format = "DNA",
#'                               lnc.format = "DNA", check.cds = TRUE,
#'                               ignore.illegal = FALSE)
#'
#' ### The input of cds.seq and lncRNA.seq can also be secondary structure
#' ### sequences when SS.features = FALSE, such as,
#' data(demp_SS.seq)
#' SS.seq <- demo_SS.seq
#' my_file_3 <- make_frequencies(cds.seq = SS.seq, lncRNA.seq = Seqs,
#'                               SS.features = FALSE, cds.format = "SS",
#'                               lnc.format = "DNA", check.cds = TRUE,
#'                               ignore.illegal = FALSE)
#' @export

make_frequencies <- function(cds.seq, mRNA.seq, lncRNA.seq, SS.features = FALSE,
                             cds.format = "DNA", lnc.format = "DNA", check.cds = TRUE,
                             ignore.illegal = TRUE) {

        if(SS.features & lnc.format != "SS") stop("Error: Secondary structure file is required for the selected mode.")

        message("+ Check the format of sequences. ", Sys.time(), "\n")

        if(lnc.format == "DNA") {
                lnc.DNA <- lncRNA.seq
        } else if(lnc.format == "SS") {
                lnc.DNA <- sapply(lncRNA.seq, get.Seq_DNA)
        } else stop("- Error: Wrong format name of lncRNA.")

        if(cds.format == "DNA") {
                cds.DNA <- cds.seq
        } else if(cds.format == "SS") {
                cds.DNA <- sapply(cds.seq, get.Seq_DNA)
        } else stop("- Error: Wrong format name of CDs.")


        if(SS.features) {
                message("+ Calculating frequencies of Secondary structure:")
                message("- Processing sequences.")
                lnc.acguD <- sapply(lncRNA.seq, get.Seq_acguD)
                pct.acguD <- sapply(mRNA.seq,   get.Seq_acguD)
                lnc.SS    <- sapply(lncRNA.seq, get.Seq_acguACGU)
                pct.SS    <- sapply(mRNA.seq,   get.Seq_acguACGU)

                message("- Calculating frequencies of acguD k = 4", "\n")

                message("  Non-coding sequences")
                lnc.acguD.freq <- get.freq(lnc.acguD, alphabet = c("a", "c", "g", "u", "D"),
                                           wordsize = 4, step = 1, freq = T, ignore.illegal = FALSE)
                message("  Completed.", "\n")

                message("  Coding sequences")
                pct.acguD.freq <- get.freq(pct.acguD, alphabet = c("a", "c", "g", "u", "D"),
                                           wordsize = 4, step = 1, freq = T, ignore.illegal = FALSE)
                message("  Completed.", "\n")

                message("- Calculating frequencies of acguACGU k = 3", "\n")

                message("  Non-coding sequences")
                lnc.acguACGU.freq <- get.freq(lnc.SS, alphabet = c("A", "C", "G", "U", "a", "c", "g", "u"),
                                              wordsize = 3, step = 1, freq = T, ignore.illegal = FALSE)
                message("  Completed.", "\n")

                message("  Coding sequences")
                pct.acguACGU.freq <- get.freq(pct.SS, alphabet = c("A", "C", "G", "U", "a", "c", "g", "u"),
                                              wordsize = 3, step = 1, freq = T, ignore.illegal = FALSE)
                message("  Completed.", "\n")
        }

        message("+ Calculating frequencies of DNA k = 6")
        if(check.cds) {
                message("- Processing CDs sequences.")
                orf.seq <- lapply(cds.DNA, function(x) {
                        orf.info <- find_ORF(x)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                cds.DNA <- orf.seq[!is.na(orf.seq)]
                message("- Calculating frequencies.")
        }

        message("\n", "  Non-coding sequences")
        lnc.DNA.freq <- get.freq(lnc.DNA, alphabet = c("a", "c", "g", "t"), wordsize = 6,
                                 step = 3, freq = F, ignore.illegal = TRUE)
        message("  Completed.", "\n")

        message("  Coding sequences")
        cds.DNA.freq <- get.freq(cds.DNA, alphabet = c("a", "c", "g", "t"), wordsize = 6,
                                 step = 3, freq = F, ignore.illegal = TRUE)
        message("  Completed.", "\n")

        message("+ Frequencies Calculation completed. ", Sys.time())

        if(SS.features) {
                Internal.freq <- list(DNA.lnc      = lnc.DNA.freq,
                                      DNA.cds      = cds.DNA.freq,
                                      acguD.lnc    = lnc.acguD.freq,
                                      acguD.pct    = pct.acguD.freq,
                                      acguACGU.lnc = lnc.acguACGU.freq,
                                      acguACGU.pct = pct.acguACGU.freq)
        } else {
                Internal.freq <- list(DNA.lnc = lnc.DNA.freq,
                                      DNA.cds = cds.DNA.freq)
        }
        Internal.freq
}

###### extract_features ######

#' Extract the Features.
#' @description This function can construct the dataset. This function is only used
#' to extract the features, please use function \code{\link{build_model}} to build
#' new models.
#'
#' @param Sequences mRNA sequences or long non-coding sequences. Can be a FASTA
#' file loaded by package "seqinr" (\code{\link[seqinr]{seqinr-package}}) or
#' secondary structure sequences (Dot-Bracket Notation) obtained from function
#' \code{\link{run_RNAfold}}. If \code{Sequences} are secondary structure
#' sequences file, parameter \code{format} should be defined as \code{"SS"}.
#'
#' @param Label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#'
#' @param SS.features Logical. If \code{SS.features = TRUE}, secondary structure
#' features will be extracted. In this case, \code{Sequences} should be secondary
#' structure sequences (Dot-Bracket Notation) obtained from function
#' \code{\link{run_RNAfold}} and parameter \code{format} should be set as \code{"SS"}.
#'
#' @param format String. Can be \code{"DNA"} or \code{"SS"}. Define the format of
#' \code{Sequences}. \code{"DNA"} for DNA sequences and \code{"SS"} for secondary
#' structure sequences. This parameter must be set as \code{"SS"} when
#' \code{SS.features = TURE}.
#'
#' @param frequencies.file String or a list obtained from function
#' \code{\link{make_frequencies}}. Input species name \code{"human"}, \code{"mouse"}
#' or \code{"wheat"} to use prebuild frequencies files. Or assign a users' own
#' frequencies file (See function \code{\link{make_frequencies}}).
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return Returns a data.frame. 11 features when \code{SS.features = FALSE} and
#' 19 features when \code{SS.features = TRUE}.
#'
#' @author Han Siyu
#' @details This function extracts the features and constructs the dataset.
#'
#' Considering that it is time consuming to obtain secondary structure sequences,
#' users can build the model only with features of sequence and EIIP
#' (\code{SS.features = FALSE}). When \code{SS.features = TRUE}, \code{Sequences}
#' should be secondary structure sequences (Dot-Bracket Notation) obtained from
#' function \code{\link{run_RNAfold}} and parameter \code{format} should be set
#' as \code{"SS"}.
#'
#' Please note that:
#'
#' Secondary structure features (\code{SS.features}) can improve the performance
#' when the species of unevaluated sequences is identical to the species of the
#' sequences that used to build the model.
#'
#' However, if users are trying to predict sequences with the model trained on
#' other species, \code{SS.features = TRUE} may lead to low accuracy.
#'
#' @section Features:
#' 1. Features based on sequence:
#'
#'    The length and coverage of the longest ORF (\code{ORF.Max.Len} and
#'    \code{ORF.Max.Cov});
#'
#'    Log-Distance.lncRNA (\code{Seq.lnc.Dist});
#'
#'    Log-Distance.protein-coding transcripts (\code{Seq.pct.Dist});
#'
#'    Distance-Ratio.sequence (\code{Seq.Dist.Ratio}).
#'
#'
#' 2. Features based on EIIP (electron-ion interaction pseudopotential) value:
#'
#'    Signal as 1/3 position (\code{Signal.Peak});
#'
#'    Signal to noise ratio (\code{SNR});
#'
#'    the minimum value of the top 10\% power spectrum (\code{Signal.Min});
#'
#'    the quantile Q1 and Q2 of the top 10\% power spectrum (\code{Singal.Q1}
#'    and \code{Signal.Q2})
#'
#'    the maximum value of the top 10\% power spectrum (\code{Signal.Max}).
#'
#'
#' 3. Features based on secondary structure sequence:
#'
#'    Log-Distance.acguD.lncRNA (\code{Dot_lnc.dist});
#'
#'    Log-Distance.acguD.protein-coding transcripts (\code{Dot_pct.dist});
#'
#'    Distance-Ratio.acguD (\code{Dot_Dist.Ratio});
#'
#'    Log-Distance.acgu-ACGU.lncRNA (\code{SS.lnc.dist});
#'
#'    Log-Distance.acgu-ACGU.protein-coding transcripts (\code{SS.pct.dist});
#'
#'    Distance-Ratio.acgu-ACGU (\code{SS.Dist.Ratio});
#'
#'    Minimum free energy (\code{MFE});
#'
#'    Percentage of Unpair-Pair (\code{UP.PCT})
#'
#' @section References:
#' HAN Siyu, LIANG Yanchun and LI Ying*.
#' LncFinder: Long Non-coding RNA Identification Tool Based on Features of
#' Sequence Intrinsic Composion, Secondary Structure and EIIP Values. (2017)
#' (\emph{Submitted})
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @seealso \code{\link{svm_tune}}, \code{\link{build_model}},
#'          \code{\link{make_frequencies}}, \code{\link{run_RNAfold}}.
#' @examples
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' ### Extract features with prebuild frequencies.file:
#' my_features <- extract_features(Seqs, Label = "Class.of.the.Sequences",
#'                                 SS.features = FALSE, format = "DNA",
#'                                 frequencies.file = "mouse",
#'                                 parallel.cores = 2)
#'
#' ### Use your own frequencies file by assign frequencies list to parameter
#' ### "frequencies.file".
#' @export

extract_features <- function(Sequences, Label = NULL, SS.features = FALSE, format = "DNA",
                             frequencies.file = "human", parallel.cores = 2) {

        if(SS.features & format != "SS") stop("Error: Secondary structure file is required for the selected mode.")

        message("+ Check the format of sequences. ", Sys.time(), "\n")
        if(format == "DNA") {
                DNA.seq <- Sequences
        } else if(format == "SS") {
                DNA.seq   <- sapply(Sequences, get.Seq_DNA)
        } else stop("- Error: Wrong format name.")

        if(class(frequencies.file) == "character") {
                if(frequencies.file == "human") {
                        Internal.data = Internal.human
                } else if(frequencies.file == "mouse") {
                        Internal.data = Internal.mouse
                } else if(frequencies.file == "wheat") {
                        Internal.data = Internal.wheat
                }
        }
         else  Internal.data = frequencies.file

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!")

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = c("Internal.data", "find_ORF"), envir = environment())

        if(SS.features) {
                message("+ Extracting secondary structure features:")
                message("- Processing sequences.")
                UP.seq    <- sapply(Sequences, get.Seq_UP)
                SS.seq    <- sapply(Sequences, get.Seq_acguACGU)
                acguD.seq <- sapply(Sequences, get.Seq_acguD)

                message("- LogDist.acguD k = 4")
                acguD    <- parallel::parSapply(cl, acguD.seq, LogDist.acguD,
                                                wordsize = 4, step = 1,
                                                hexamer.lnc = Internal.data$acguD.lnc,
                                                hexamer.cds = Internal.data$acguD.pct)

                message("- LogDist.acguACGU k = 3")
                acguACGU <- parallel::parSapply(cl, SS.seq, LogDist.acguACGU,
                                                wordsize = 3, step = 1,
                                                hexamer.lnc = Internal.data$acguACGU.lnc,
                                                hexamer.cds = Internal.data$acguACGU.pct)
                message("- MFE")
                MFE      <- as.numeric(as.character(Sequences[3,]))

                message("- Percentage of Unpair-Pair", "\n")
                UP.percentage <- parallel::parSapply(cl, UP.seq, seqinr::count,
                                                     wordsize = 2, freq = TRUE,
                                                     alphabet = c("U", "P"))
                SS.group <- data.frame(t(acguD), t(acguACGU), MFE = MFE, UP.PCT = UP.percentage[3,])
        }

        message("+ Extracting sequence and EIIP features:")
        message("- ORF and LogDist.Seq k = 3")
        Seq.ORF  <- parallel::parSapply(cl, DNA.seq, LogDist.DNA, wordsize = 6, step = 3,
                                        hexamer.lnc = Internal.data$DNA.lnc,
                                        hexamer.cds = Internal.data$DNA.cds)

        message("- EIIP", "\n")
        EIIP.DFT <- parallel::parSapply(cl, DNA.seq, EIIP.DFT)
        parallel::stopCluster(cl)
        features <- data.frame(t(Seq.ORF), t(EIIP.DFT))

        message("+ Feature extraction completed. ", Sys.time())
        if(SS.features) features <- cbind(features, SS.group)
        if(!is.null(Label)) features <- cbind(Label, features)
        features
}

###### svm_tune ######

#' Parameter Tuning of SVM
#' @description This function conduct the parameter tuning of SVM. Parameters
#' gamma and cost can be tuned using grid search.
#'
#' @param dataset The dataset obtained from function \code{\link{extract_features}}.
#'
#' @param folds.num Integer. Specify the number of folds for cross-validation.
#' (Default: \code{10})
#'
#' @param seed Integer. Used to set the seed for cross-validation. (Default: \code{1})
#'
#' @param gamma.range The range of gamma. (Default: \code{2 ^ seq(-5, 0, 1)})
#'
#' @param cost.range The range of cost. (Default: \code{c(1, 4, 8, 16, 24, 32)})
#'
#' @param return.model Logical. If \code{TRUE}, the function will return the best
#' model trained on the full dataset. If \code{FALSE}, this function will return
#' the optimal parameters.
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}, users can set as \code{-1} to run
#' this function with all cores. If the number of \code{parallel.cores} is more
#' than the \code{folds.num} (number of the folds for cross-validation), the
#' number of \code{parallel.cores} will be set as \code{folds.num} automatically.
#'
#' @return Returns the optimal parameters when \code{return.model = FALSE}.
#'
#' Or returns the best model when \code{return.model = TRUE}.
#' @author Han Siyu
#'
#' @details During the model tuning, the performance of each combination of
#' parameters will output. Sensitivity, Specificity, Accuracy, F-Measure and Kappa
#' Value are used to evaluate the performances. The best gamma and cost (or best
#' model) are selected based on Accuracy.
#'
#' For the details of parameter gamma and cost, please refer to function
#' \code{\link[e1071]{svm}} of package "e1071".
#'
#' For the details of metrics, please refer to function
#' \code{\link[caret]{confusionMatrix}} of package "caret".
#'
#' @importFrom caret createFolds
#' @importFrom caret confusionMatrix
#' @importFrom e1071 svm
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @seealso extract_features
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' positive_data <- extract_features(Seqs[1:5], Label = "NonCoding",
#'                                   SS.features = FALSE, format = "DNA",
#'                                   frequencies.file = "human",
#'                                   parallel.cores = 2)
#'
#' negative_data <- extract_features(Seqs[6:10], Label = "Coding",
#'                                   SS.features = FALSE, format = "DNA",
#'                                   frequencies.file = "human",
#'                                   parallel.cores = 2)
#'
#' my_dataset <- rbind(positive_data, negative_data)
#'
#' ### Or use our data "demo_dataset"
#' data(demo_dataset)
#' my_dataset <- demo_dataset
#'
#' optimal_parameter <- svm_tune(my_dataset, folds.num = 2, seed = 1,
#'                               gamma.range = (2 ^ seq(-5, 0, 2)),
#'                               cost.range = c(1, 8, 16),
#'                               return.model = FALSE, parallel.core = 2)
#'
#' ### Users can set return.model = TRUE to return the best model.
#' }
#' @export

svm_tune <- function(dataset, folds.num = 10, seed = 1,
                     gamma.range = (2 ^ seq(-5, 0, 1)),
                     cost.range = c(1, 4, 8, 16, 24, 32),
                     return.model = TRUE, parallel.cores = 2){

        set.seed(seed)
        folds <- caret::createFolds(dataset$Label, k = folds.num, returnTrain = TRUE)

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        parallel.cores <- ifelse(parallel.cores > folds.num, folds.num, parallel.cores)

        if(parallel.cores == 2 & folds.num > parallel.cores) {
                message("Users can try to set parallel.cores = -1 to use all cores!")
        }

        cl <- parallel::makeCluster(parallel.cores)

        res <- c()
        message("+ SVM.tune Processing.")
        for(g in gamma.range){
                for(C in cost.range){
                        parallel::clusterExport(cl, varlist = c("g", "C", "dataset"),
                                                envir = environment())
                        message("- gamma = ", g, ", Cost = ", C)
                        perf.res    <- parallel::parSapply(cl, folds, function(x, gamma = g, cost = C) {
                                subset <- dataset[x, ]
                                test.set  <- dataset[-x, ]
                                svm.mod <- e1071::svm(Label ~ ., data = subset, scale = TRUE, probability = TRUE,
                                               kernel = "radial", gamma = gamma, cost = cost)

                                res <- stats::predict(svm.mod, test.set, probability = TRUE)

                                prob <- data.frame(attr(res, "probabilities"), Label = test.set$Label)

                                confusion.res <- caret::confusionMatrix(data.frame(res)$res, test.set$Label,
                                                                        positive = "NonCoding", mode = "everything")
                                performance.res <- data.frame(Sensitivity = confusion.res$byClass[1],
                                                              Specificity = confusion.res$byClass[2],
                                                              Accuracy    = confusion.res$overall[1],
                                                              F.Measure   = confusion.res$byClass[7],
                                                              Kappa       = confusion.res$overall[2])

                                performance.res
                        })
                        Ave.res     <- apply(perf.res, 1, as.numeric)
                        Ave.res     <- as.data.frame(t(Ave.res))
                        Ave.res$Ave.Res <- rowMeans(Ave.res)
                        print(t(Ave.res[ncol(Ave.res)]))
                        res <- c(res, list(Ave.res))
                        names(res)[length(res)] <- paste("g = ", g, "; C = ", C, sep = "")
                }
        }
        parallel::stopCluster(cl)
        message("\n", "Best Result:")
        Accuracy <- sapply(res, function(x) x[3,ncol(Ave.res)])
        best.res <- res[Accuracy == max(Accuracy)][1]
        message("$ ", names(best.res))
        print(t(best.res[[1]][ncol(Ave.res)]))
        best.gamma <- as.numeric(strsplit(strsplit(names(best.res), "; ")[[1]][[1]], " = ")[[1]][[2]])
        best.cost  <- as.numeric(strsplit(strsplit(names(best.res), "; ")[[1]][[2]], " = ")[[1]][[2]])
        best.parameters <- c(best.gamma = best.gamma, best.cost = best.cost)
        if(return.model) {
                message("\n", "+ Training the Model on the Full Dataset.")
                svm.mod  <- e1071::svm(Label ~ ., data = dataset, scale = TRUE,
                                       probability = TRUE, kernel = "radial",
                                       gamma = best.parameters[[1]],
                                       cost  = best.parameters[[2]])
                message("\n", "+ SVM.tune completed.")
                return(svm.mod)
        } else {
                message("\n", "+ SVM.tune completed.")
                return(best.parameters)
        }
}

###### build_model ######

#' Build Users' Own Model
#' @description This function is used to build new models with users' own data.
#'
#' @param lncRNA.seq Long non-coding sequences. Can be a FASTA file loaded by
#' package "seqinr" (\code{\link[seqinr]{seqinr-package}}) or secondary structure
#' sequences file (Dot-Bracket Notation) obtained from function
#' \code{\link{run_RNAfold}}. If \code{lncRNA.seq} is secondary structure
#' sequences file, parameter \code{lncRNA.format} should be defined as \code{"SS"}.
#'
#' @param mRNA.seq mRNA sequences. FASTA file loaded by package "seqinr" or
#' secondary structure sequences (Dot-Bracket Notation) obtained from function
#' \code{\link{run_RNAfold}}. If \code{mRNA.seq} is secondary structure sequences
#' file, parameter \code{mRNA.format} should be defined as \code{"SS"}.
#'
#' @param frequencies.file String or a list obtained from function
#' \code{\link{make_frequencies}}. Input species name \code{"human"},
#' \code{"mouse"} or \code{"wheat"} to use prebuild frequencies files. Or assign
#' a users' own frequencies file (Please refer to function
#' \code{\link{make_frequencies}} for more information).
#'
#' @param SS.features Logical. If \code{SS.features = TRUE}, secondary structure
#' features will be used to build the model. In this case, \code{lncRNA.seq} and
#' \code{mRNA.seq} should be secondary structure sequences (Dot-Bracket Notation)
#' obtained from function \code{\link{run_RNAfold}} and parameter
#' \code{lncRNA.format} and \code{mRNA.format} should be set as \code{"SS"}.
#'
#' @param lncRNA.format String. Define the format of \code{lncRNA.seq}. \code{"DNA"}
#' for DNA sequences and \code{"SS"} for secondary structure sequences. Only when
#' both \code{mRNA.format} and \code{lncRNA.format} are aet as \code{"SS"}, can
#' the model with secondary structure features be built (\code{SS.features = TRUE}).
#'
#' @param mRNA.format String. Define the format of \code{mRNA.seq}. Can be
#' \code{"DNA"} or \code{"SS"}. \code{"DNA"} for DNA sequences and \code{"SS"}
#' for secondary structure sequences. When this parameter is defined as \code{"DNA"},
#' only the model without secondary structure features can be built. In this case,
#' parameter \code{SS.features} should be set as \code{FALSE}.
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}, users can set as \code{-1} to run
#' this function with all cores. During the process of svm tuning, if the number
#' of \code{parallel.cores} is more than the \code{folds.num} (number of the folds
#' for cross-validation), the number of \code{parallel.cores} will be set as
#' \code{folds.num} automatically.
#'
#' @param folds.num Integer. Specify the number of folds for cross-validation.
#' (Default: \code{10})
#'
#' @param seed Integer. Used to set the seed for cross-validation. (Default: \code{1})
#'
#' @param gamma.range The range of gamma. (Default: \code{2 ^ seq(-5, 0, 1)})
#'
#' @param cost.range The range of cost. (Default: \code{c(1, 4, 8, 16, 24, 32)})
#'
#' @return Returns a svm model.
#'
#' @author Han Siyu
#'
#' @details This function is used to build a new model with users' own sequences.
#' Users can use function \code{\link{lnc_finder}} to predict the sequences with
#' new models.
#'
#' For the details of \code{frequencies.file}, please refer to function
#' \code{\link{make_frequencies}}.
#'
#' For the details of the features, please refer to function
#' \code{\link{extract_features}}.
#'
#' For the details of svm tuning, please refer to function \code{\link{svm_tune}}.
#'
#' @section References:
#' HAN Siyu, LIANG Yanchun and LI Ying*.
#' LncFinder: Long Non-coding RNA Identification Tool Based on Features of
#' Sequence Intrinsic Composion, Secondary Structure and EIIP Values. (2017)
#' (\emph{Submitted})
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom caret createFolds
#' @importFrom caret confusionMatrix
#' @importFrom e1071 svm
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @seealso \code{\link{make_frequencies}}, \code{\link{lnc_finder}},
#'          \code{\link{extract_features}}, \code{\link{svm_tune}},
#'          \code{\link[e1071]{svm}}.
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' ### Build the model with prebuild frequencies.file:
#' my_model <- build_model(lncRNA.seq = Seqs[1:5], mRNA.seq = Seqs[6:10],
#'                         frequencies.file = "human", SS.features = FALSE,
#'                         lncRNA.format = "DNA", mRNA.format = "DNA",
#'                         parallel.cores = 2, folds.num = 2, seed = 1,
#'                         gamma.range = (2 ^ seq(-5, -1, 2)),
#'                         cost.range = c(2, 6, 12, 20))
#'
#' ### Users can use default values of gamma.range and cost.range to find the
#' best parameters.
#' ### Use your own frequencies file by assigning frequencies list to parameter
#' ### "frequencies.file".
#' }
#' @export

build_model <- function(lncRNA.seq, mRNA.seq, frequencies.file, SS.features = FALSE,
                        lncRNA.format = "DNA", mRNA.format = "DNA", parallel.cores = 2,
                        folds.num = 10, seed = 1, gamma.range = (2 ^ seq(-5, 0, 1)),
                        cost.range = c(1, 4, 8, 16, 24, 32)) {

        message("+ Extract features of non-coding sequences.")
        lnc.data <- extract_features(lncRNA.seq, Label = "NonCoding", SS.features = SS.features,
                                     format = lncRNA.format, frequencies.file = frequencies.file,
                                     parallel.cores = parallel.cores)

        message("+ Extract features of protein-coding sequences.")
        pct.data <- extract_features(mRNA.seq, Label = "Coding", SS.features = SS.features,
                                     format = mRNA.format, frequencies.file = frequencies.file,
                                     parallel.cores = parallel.cores)

        dataset  <- rbind(lnc.data, pct.data)

        best.parameters <- svm_tune(dataset, gamma.range = gamma.range, cost.range = cost.range,
                                    seed = seed, folds.num = folds.num, return.model = FALSE,
                                    parallel.cores = parallel.cores)

        message("+ Build the model with the whole dataset.")
        svm.mod  <- e1071::svm(Label ~ ., data = dataset, scale = TRUE, probability = TRUE, kernel = "radial",
                        gamma = best.parameters[[1]], cost = best.parameters[[2]])

        message("+ Process completed.")
        svm.mod
}

###### lnc_finder ######

#' Long Non-coding RNA Identification
#' @description This function is used to predict sequences are non-coding transcripts
#' or protein-coding transcripts.
#'
#' @param Sequences Unevaluated sequences. Can be a FASTA file loaded by package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}) or secondary structure sequences
#' (Dot-Bracket Notation) obtained from function \code{\link{run_RNAfold}}. If
#' \code{Sequences} is secondary structure sequences file, parameter \code{format}
#' should be defined as \code{"SS"}.
#'
#' @param SS.features Logical. If \code{SS.features = TRUE}, secondary structure
#' features will be used. In this case, \code{Sequences} should be secondary
#' structure sequences (Dot-Bracket Notation) obtained from function
#' \code{\link{run_RNAfold}} and parameter \code{format} should be set as \code{"SS"}.
#'
#' @param format String. Define the format of the \code{Sequences}. Can be
#' \code{"DNA"} or \code{"SS"}. \code{"DNA"} for DNA sequences and \code{"SS"}
#' for secondary structure sequences. This parameter must be set as \code{"SS"}
#' when \code{SS.features = TURE}.
#'
#' @param frequencies.file String or a list obtained from function
#' \code{\link{make_frequencies}}. Input species name \code{"human"}, \code{"mouse"}
#' or \code{"wheat"} to use prebuild frequencies files. Or assign a users' own
#' frequencies file (See function \code{\link{make_frequencies}}).
#'
#' @param svm.model String or a svm model obtained from funtion \code{\link{build_model}}
#' or \code{\link{svm_tune}}. Input species name \code{"human"}, \code{"mouse"}
#' or \code{"wheat"} to use prebuild models. Or assign a users' own model (See
#' function \code{\link{build_model}}).
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return Returns a data.frame. Including the results of prediction (\code{Pred});
#' coding potential (\code{Coding.Potential}) and the features. For the details
#' of the features, please refer to function \code{\link{extract_features}}.
#'
#' @author Han Siyu
#'
#' @details Considering that it is time consuming to obtain secondary structure
#' sequences, users can input nucleotide sequences and predict these sequences
#' without secondary structure features (Set \code{SS.features} as \code{FALSE}).
#'
#' Please note that:
#'
#' \code{SS.features} can improve the performance when the species of unevaluated
#' sequences is identical to the species of the sequences that used to build the
#' model.
#'
#' However, if users are trying to predict sequences with the model trained on
#' other species, \code{SS.features} may lead to low accuracy.
#'
#' For the details of \code{frequencies.file}, please refer to function
#' \code{\link{make_frequencies}}.
#'
#' For the details of the features, please refer to function
#' \code{\link{extract_features}}.
#'
#' @section References:
#' HAN Siyu, LIANG Yanchun and LI Ying*.
#' LncFinder: Long Non-coding RNA Identification Tool Based on Features of
#' Sequence Intrinsic Composion, Secondary Structure and EIIP Values. (2017)
#' (\emph{Submitted})
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @seealso \code{\link{build_model}}, \code{\link{make_frequencies}},
#'          \code{\link{extract_features}}, \code{\link{run_RNAfold}}.
#' @examples
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' ### Input one sequence:
#' OneSeq <- Seqs[1]
#' result_1 <- lnc_finder(OneSeq, SS.features = FALSE, format = "DNA",
#'                        frequencies.file = "human", svm.model = "human",
#'                        parallel.cores = 2)
#'
#' \dontrun{
#' ### Or several sequences:
#' data(demo_SS.seq)
#' Seqs <- demo_SS.seq
#' result_2 <- lnc_finder(Seqs, SS.features = TRUE, format = "SS",
#'                        frequencies.file = "mouse", svm.model = "mouse",
#'                        parallel.cores = 2)
#'
#' ### A complete work flow:
#' ### Calculate second structure on Windows OS,
#' RNAfold.path <- '"E:/Program Files/ViennaRNA/RNAfold.exe"'
#' SS.seq <- run_RNAfold(Seqs, RNAfold.path = RNAfold.path, parallel.cores = 2)
#'
#' ### Predict the sequences with secondary structure features,
#' result_2 <- lnc_finder(SS.seq, SS.features = TRUE, format = "SS",
#'                        frequencies.file = "mouse", svm.model = "mouse",
#'                        parallel.cores = 2)
#'
#' ### Predict sequences with your own model by assigning a new svm.model and
#' ### frequencies.file to parameters "svm.model" and "frequencies.file"
#' }
#' @export

lnc_finder <- function(Sequences, SS.features = FALSE, format = "DNA", frequencies.file = "human",
                       svm.model = "human", parallel.cores = 2) {

        if(class(svm.model) == "character") {
                if(SS.features) {
                        if(svm.model == "human") {
                                svm.mod <- human.mod
                        } else if(svm.model == "mouse") {
                                svm.mod <- mouse.mod
                        } else if(svm.model == "wheat") {
                                svm.mod <- wheat.mod
                        } else stop("Wrong svm.model name.")
                } else {
                        if(svm.model == "human") {
                                svm.mod <- human.mod.no_ss
                        } else if(svm.model == "mouse") {
                                svm.mod <- mouse.mod.no_ss
                        } else if(svm.model == "wheat") {
                                svm.mod <- wheat.mod.no_ss
                        } else stop("Wrong svm.model name.")
                }
        } else svm.mod <- svm.model

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)

        seq.data <- extract_features(Sequences, SS.features = SS.features, format = format,
                                     frequencies.file = frequencies.file, parallel.cores = parallel.cores)

        pred.res <- stats::predict(svm.mod, seq.data, probability = TRUE)
        results  <- data.frame(Pred = pred.res)
        prob.res <- data.frame(Result = results, Coding.Potential = (attr(pred.res, "probabilities")[,2]))
        output   <- cbind(prob.res, seq.data)
        output
}
