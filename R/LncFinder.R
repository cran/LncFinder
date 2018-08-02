################  Internal Functions  ################

Internal.convertProb <- function(value, base, param, prob, weight) {
        for(idx in 1:10) {
                if(value >= param[idx]) {
                        probRes <- prob[[base]][[idx]] * weight[[base]]
                        break
                }
        }
        probRes
}

Internal.FickettScore <- function(oneSeq, .positionProb, .positionWeight,
                                  .positionParam, .contentProb,
                                  .contentWeight, .contentParam) {
        countSeq <- seqinr::count(oneSeq, wordsize = 1, freq = TRUE)
        lenSeq <- length(oneSeq)
        phase1 <- oneSeq[seq(1, lenSeq, 3)]
        phase2 <- oneSeq[seq(2, lenSeq, 3)]
        phase3 <- oneSeq[seq(3, lenSeq, 3)]

        position.a <- max(sum(phase1 == "a"), sum(phase2 == "a"), sum(phase3 == "a")) / (min(sum(phase1 == "a"), sum(phase2 == "a"), sum(phase3 == "a")) + 1)
        position.c <- max(sum(phase1 == "c"), sum(phase2 == "c"), sum(phase3 == "c")) / (min(sum(phase1 == "c"), sum(phase2 == "c"), sum(phase3 == "c")) + 1)
        position.g <- max(sum(phase1 == "g"), sum(phase2 == "g"), sum(phase3 == "g")) / (min(sum(phase1 == "g"), sum(phase2 == "g"), sum(phase3 == "g")) + 1)
        position.t <- max(sum(phase1 == "t"), sum(phase2 == "t"), sum(phase3 == "t")) / (min(sum(phase1 == "t"), sum(phase2 == "t"), sum(phase3 == "t")) + 1)

        contentProb.a <- Internal.convertProb(value = countSeq["a"], base = "a", param = .contentParam, prob = .contentProb, weight = .contentWeight)
        contentProb.c <- Internal.convertProb(value = countSeq["c"], base = "c", param = .contentParam, prob = .contentProb, weight = .contentWeight)
        contentProb.g <- Internal.convertProb(value = countSeq["g"], base = "g", param = .contentParam, prob = .contentProb, weight = .contentWeight)
        contentProb.t <- Internal.convertProb(value = countSeq["t"], base = "t", param = .contentParam, prob = .contentProb, weight = .contentWeight)

        positionProb.a <- Internal.convertProb(value = position.a, base = "a", param = .positionParam, prob = .positionProb, weight = .positionWeight)
        positionProb.c <- Internal.convertProb(value = position.c, base = "c", param = .positionParam, prob = .positionProb, weight = .positionWeight)
        positionProb.g <- Internal.convertProb(value = position.g, base = "g", param = .positionParam, prob = .positionProb, weight = .positionWeight)
        positionProb.t <- Internal.convertProb(value = position.t, base = "t", param = .positionParam, prob = .positionProb, weight = .positionWeight)

        FickettScore <- sum(contentProb.a, contentProb.c, contentProb.g, contentProb.t,
                            positionProb.a, positionProb.c, positionProb.g, positionProb.t)
        FickettScore
}

Internal.EucDistance <- function(OneSeq, referFreq, k, step, alphabet){

        freq.val <- seqinr::count(OneSeq, wordsize = k, by = step, alphabet = alphabet, freq = TRUE)

        EucDist.LNC <- sqrt(sum((freq.val - referFreq$ref.lnc) ^ 2))
        EucDist.PCT <- sqrt(sum((freq.val - referFreq$ref.cds) ^ 2))
        EucDist.Ratio <- EucDist.LNC / EucDist.PCT

        distance <- c(EucDist.LNC = EucDist.LNC, EucDist.PCT = EucDist.PCT, EucDist.Ratio = EucDist.Ratio)

        distance
}

Internal.LogDistance <- function(OneSeq, referFreq, k, step, alphabet){

        count.num <- seqinr::count(OneSeq, wordsize = k, by = step, alphabet = alphabet, freq = FALSE)

        freq.val  <- count.num / sum(count.num)
        lnc.ratio <- freq.val[referFreq$ref.lnc != 0] / referFreq$ref.lnc[referFreq$ref.lnc != 0]
        pct.ratio <- freq.val[referFreq$ref.cds != 0] / referFreq$ref.cds[referFreq$ref.cds != 0]
        LogDist.LNC <- (sum(log(lnc.ratio[lnc.ratio != 0])) - length(which(lnc.ratio == 0)) + length(which(referFreq$ref.lnc == 0))) / sum(count.num)
        LogDist.PCT <- (sum(log(pct.ratio[pct.ratio != 0])) - length(which(pct.ratio == 0)) + length(which(referFreq$ref.cds == 0))) / sum(count.num)
        LogDist.Ratio <- LogDist.LNC / LogDist.PCT

        distance <- c(LogDist.LNC = LogDist.LNC, LogDist.PCT = LogDist.PCT, LogDist.Ratio = LogDist.Ratio)

        distance
}

Internal.hexamerScore <- function(OneSeq, referFreq, k, step, alphabet){

        count.num <- seqinr::count(OneSeq, wordsize = k, by = step, alphabet = alphabet, freq = FALSE)

        seq.sum  <- 0
        for(i in names(count.num)[count.num != 0]){
                if(referFreq$ref.cds[[i]] != 0 & referFreq$ref.cds[[i]] != 0) {
                        seq.sum <- seq.sum + (count.num[[i]] * log(referFreq$ref.cds[[i]] / referFreq$ref.cds[[i]]))
                } else {
                        if(referFreq$ref.cds[[i]] == 0) {
                                seq.sum <- seq.sum - count.num[[i]]
                        }
                        if(referFreq$ref.lnc[[i]] == 0) {
                                seq.sum <- seq.sum + count.num[[i]]
                        }
                }
        }
        hex.score <- seq.sum / sum(count.num)

        hex.score
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
        ordered.pow  <- sort(seq.spectrum, decreasing = TRUE)
        pow.percent  <- ordered.pow[1:round(length(ordered.pow) * 0.1)]
        Quantile.res <- stats::quantile(pow.percent)[-4]
        names(Quantile.res) <- c("Signal.Min", "Signal.Q1", "Signal.Q2", "Signal.Max")
        EIIP.res <- c(SNR.res, Quantile.res)

        return(EIIP.res)
}

get_EIIP <- function(OneSeq, percent.length = 0.1, probs = seq(0, 1, 0.25)) {

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
        SNR.res <- c(Signal.Peak = DFT.max, Average.Power = total.power, SNR = (DFT.max / total.power))

        seq.spectrum <- seq.spectrum[-c(1, power.length)]

        ordered.pow  <- sort(seq.spectrum, decreasing = TRUE)
        pow.percent  <- ordered.pow[1:round(length(ordered.pow) * percent.length)]
        Quantile.res <- stats::quantile(pow.percent, probs = probs)

        EIIP.res <- c(SNR.res, Quantile.res)

        return(EIIP.res)
}

fold.res <- function(x, dataset, positive.class, ...) {
        subset <- dataset[x, ]
        test.set  <- dataset[-x, ]

        svm.mod <- e1071::svm(label ~ ., data = subset, probability = TRUE, ...)
        print(svm.mod)

        res <- stats::predict(svm.mod, test.set, probability = TRUE)

        confusion.res <- caret::confusionMatrix(data.frame(res)[,1], test.set$label,
                                                positive = positive.class, mode = "everything")
        performance.res <- data.frame(Sensitivity = confusion.res$byClass[1],
                                      Specificity = confusion.res$byClass[2],
                                      Accuracy    = confusion.res$overall[1],
                                      F.Measure   = confusion.res$byClass[7],
                                      Kappa       = confusion.res$overall[2])
        performance.res
}

get.svm.res <- function(x, index, dataset = dataset){

        train.set <- dataset[x,c("label", index)]
        test.set  <- dataset[-x,]

        svm.mod   <- e1071::svm(label ~ ., data = train.set, scale = TRUE, probability = TRUE, kernel = "radial")
        results   <- stats::predict(svm.mod, test.set, probability = TRUE)
        conf.res  <- caret::confusionMatrix(data.frame(results)[,1], test.set$label,
                                            positive = "lnc", mode = "everything")
        perf.res  <- data.frame(Sensitivity = conf.res$byClass[1],
                                Specificity = conf.res$byClass[2],
                                Accuracy    = conf.res$overall[1],
                                F.Measure   = conf.res$byClass[7],
                                Kappa       = conf.res$overall[2])

        w.value <- t(svm.mod$coefs) %*% svm.mod$SV
        weights <- w.value * w.value
        weights <- data.frame(weight = t(weights))

        output  <- list(perf.res, weights)
        output
}

# svm.rfe <- function(dataset, variable.group, save.rank = FALSE, file.path = NULL){
#         message("+ Initialization.  ",  Sys.time())
#         message("")
#
#         if(is.null(file.path)) file.path <- getwd()
#
#         if(save.rank) setwd(file.path)
#         variable.group <- c(variable.group, 0)
#         variable.group <- unique(variable.group)
#         variable.group <- variable.group[which(variable.group < (ncol(dataset) - 1) & variable.group >= 0)]
#         variable.group <- sort(variable.group, decreasing = TRUE)
#         feature.num    <- ncol(dataset)
#
#         set.seed(1)
#         folds <- caret::createFolds(dataset$label, k = 10, returnTrain = TRUE)
#
#         new.fea <- c()
#
#         perf.all <- sapply(variable.group, function(variable.number){
#
#                 varindex <- if(variable.number == max(variable.group)) names(dataset)[-1] else new.fea
#                 message("+ RFE Processing... Variables: ", length(varindex))
#                 cl <- parallel::makeCluster(10)
#                 parallel::clusterExport(cl, varlist = c("dataset", "get.svm.res", "varindex"), envir = environment())
#                 message("- Training Model...  (",  Sys.time(), ")")
#                 svm.res <- parallel::parLapply(cl, folds, get.svm.res, index = varindex)
#                 parallel::stopCluster(cl)
#                 message("- Completed.         (",  Sys.time(), ")")
#                 message("- Metric:")
#
#                 perf    <- sapply(svm.res, rbind)
#                 perf    <- sapply(perf[1,], cbind)
#                 perf.convert <- apply(perf, 1, as.numeric)
#                 print(colMeans(perf.convert))
#
#                 weights <- svm.res[[1]][[2]]
#                 for(weight.index in 2:10) {
#                         weight  <- svm.res[[weight.index]][[2]]
#                         weights <- weights + weight
#                 }
#                 weights    <- (weights - min(weights)) / (max(weights) - min(weights)) * 100
#                 sort.index <- sort(weights$weight, decreasing = T, index.return = TRUE)$ix
#                 weights    <- weights[sort.index, , drop = F]
#                 if(nrow(weights) > 1 & save.rank) {
#                         save(weights, file = paste("ranking.", nrow(weights), ".RData", sep = ""))
#                 }
#                 if(variable.number != 0) {
#                         new.fea    <<- rownames(weights)[1:variable.number]
#                         message("- Remained Features:")
#                         cat(new.fea, fill = TRUE, labels = " ", sep = "  ")
#                         message("- Omitted Features:")
#                         cat(rownames(weights)[-(1:variable.number)], fill = TRUE, labels = " ", sep = "  ")
#                         message("- New Dataset Loaded. Variables: ", variable.number)
#                         message("- Remained Iterations: ", length(variable.group) - which(variable.group == variable.number))
#                         message("")
#                 }
#                 perf
#         })
#
#         ### Assign Column Names
#         var.colnames <- variable.group[which(variable.group != 0)]
#         var.colnames <- c(feature.num - 1, var.colnames)
#         colnames(perf.all) <- var.colnames
#
#         ### Assign Row Names
#         folds.names    <- sapply(c(1:10), function(x) paste("Fold", x, sep = ""))
#         rfe.rownames   <- c()
#         for(folds.index in folds.names){
#                 for (metric in c("Sensitivity", "Specificity", "Accuracy", "F.Measure", "Kappa")) {
#                         temp.name    <- paste(metric, folds.index, sep = ".")
#                         rfe.rownames <- c(rfe.rownames, temp.name)
#                 }
#         }
#         row.names(perf.all) <- rfe.rownames
#
#         ### Obtain the Average Performance
#         Ave <- apply(perf.all, 2, function(col){
#                 Ave.res <- data.frame(Sensitivity = mean(unlist(col[seq(1, 50, 5)])),
#                                       Specificity = mean(unlist(col[seq(2, 50, 5)])),
#                                       Accuracy    = mean(unlist(col[seq(3, 50, 5)])),
#                                       F.Measure   = mean(unlist(col[seq(4, 50, 5)])),
#                                       Kappa       = mean(unlist(col[seq(5, 50, 5)])))
#
#         })
#         Ave.res <- sapply(Ave, rbind)
#         Ave.res <- t(Ave.res)
#         Ave.res <- as.data.frame(apply(Ave.res, 2, unlist))
#         # Ave.res <- as.data.frame(Ave.res[nrow(Ave.res):1,])
#
#         perf.all.final <- list(raw.perf = perf.all, ave.perf = Ave.res)
#
#         message("")
#         message("+ RFE Process Completed. ", Sys.time())
#         perf.all.final
# }

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

find_ORF <- function(OneSeq, max.only = TRUE) {
        OneSeq <- unlist(seqinr::getSequence(OneSeq, as.string = TRUE))
        OneSeq <- gsub("\n", "", OneSeq)
        start_pos <- unlist(gregexpr("atg", OneSeq, ignore.case = TRUE))
        if(sum(start_pos) == -1) {
                if(max.only) {
                        ORF_info <- list(ORF.Max.Len = 0, ORF.Max.Cov = 0, ORF.Max.Seq = NA)
                } else {
                        ORF_info <- data.frame(ORF.Seq = NA, ORF.Start = 0, ORF.Stop = 0, ORF.Len = 0)
                }

        } else {
                stop_pos <- sapply(c("taa", "tag", "tga"), function(x){
                        pos <- unlist(gregexpr(x, OneSeq, ignore.case = TRUE))
                })
                stop_pos <- sort(unlist(stop_pos, use.names = FALSE))
                seq_length <- nchar(OneSeq)
                if(all(stop_pos == -1)) {
                        if(max.only) {
                                orf_start <- min(start_pos)
                                orf_stop  <- seq_length - ((seq_length - orf_start + 1) %% 3) - 2
                                max_len   <- orf_stop - orf_start + 3
                                max_cov   <- max_len / seq_length
                                orf_seq   <- substr(OneSeq, start = orf_start, stop = orf_stop + 2)
                                ORF_info  <- list(ORF.Max.Len = max_len, ORF.Max.Cov = max_cov, ORF.Max.Seq = orf_seq)
                        } else {
                                for(i in seq_along(start_pos)) {
                                        start.val  <- start_pos[i]
                                        stop.val   <- seq_length - ((seq_length - start.val + 1) %% 3) - 2
                                        length.val <- stop.val - start.val + 3
                                        cover.val  <- length.val / seq_length
                                        ORF.seq    <- substr(OneSeq, start = start.val, stop = stop.val + 2)
                                        output.tmp <- data.frame(ORF.Seq = ORF.seq, ORF.Start = start.val, ORF.Stop = stop.val,
                                                                 ORF.Len = length.val, ORF.Cov = cover.val, stringsAsFactors = FALSE)
                                        if(i == 1) {
                                                ORF_info <- output.tmp
                                        } else {
                                                ORF_info <- rbind(ORF_info, output.tmp)
                                        }

                                }
                        }

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
                        if(max.only) {
                                max_len   <- max(orf_lengths)
                                max_cov   <- max_len / seq_length
                                max_index <- which(orf_lengths == max_len)
                                orf_start <- orf_starts[max_index]
                                orf_stop  <- orf_stops[max_index]
                                orf_seq   <- substr(OneSeq, start = orf_start, stop = orf_stop + 2)
                                ORF_info  <- list(ORF.Max.Len = max_len, ORF.Max.Cov = max_cov, ORF.Max.Seq = orf_seq)
                        } else {
                                for(i in seq_along(orf_starts)) {
                                        start.val  <- orf_starts[i]
                                        stop.val   <- orf_stops[i] + 2
                                        length.val <- orf_lengths[i]
                                        cover.val  <- length.val / seq_length
                                        ORF.seq    <- substr(OneSeq, start = start.val, stop = stop.val)
                                        output.tmp <- data.frame(ORF.Seq = ORF.seq, ORF.Start = start.val, ORF.Stop = stop.val,
                                                                 ORF.Len = length.val, ORF.Cov = cover.val, stringsAsFactors = FALSE)
                                        if(i == 1) {
                                                ORF_info <- output.tmp
                                        } else {
                                                ORF_info <- rbind(ORF_info, output.tmp)
                                        }

                                }
                        }
                }

        }
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

secondary_seq <- function(OneSeq, info, RNAfold.path){
        seq.string <- unlist(seqinr::getSequence(OneSeq, TRUE))

        # assign("index", index + 1, inherits = TRUE)
        index <- get("index")
        showMessage <- paste(index, "of", info, "length:", nchar(seq.string), "nt", "\n")
        cat(showMessage)

        RNAfold.command <- paste(RNAfold.path, "--noPS")

        seq.ss <- system(RNAfold.command, intern = TRUE, input = seq.string)
        index <<- index + 1

        seq.ss[3] <- as.numeric(substr(seq.ss[2], nchar(seq.string) + 3, nchar(seq.ss[2]) - 1))
        seq.ss[2] <- substr(seq.ss[2], 1, nchar(seq.string))
        seq.ss
}

# CPAT.evaluate <- function(path) {
#
#         dataset <- read.table(file = path, sep ="\t", header = T)
#
#         set.seed(1)
#         folds <- caret::createFolds(dataset$label, k = 10, returnTrain = TRUE)
#
#         threshold <- sapply(folds, function(x) {
#                 subset <- dataset[x, ]
#                 test.set  <- dataset[-x, ]
#                 glm.mod <- glm(label ~ mRNA + ORF + Fickett + Hexamer, data = subset, family = binomial(link="logit"))
#                 res <- predict(glm.mod, test.set, type = "response")
#                 res <- data.frame(Prob = res, label = test.set$label)
#                 threshold <- pROC::coords(pROC::roc(res$label, res$Prob), x = "best", ret='threshold')[[1]]
#         })
#         best.cutoff <- round(sum(threshold) / 10, digits = 4)
#
#         performance <- sapply(folds, function(x) {
#                 subset <- dataset[x, ]
#                 test.set  <- dataset[-x, ]
#                 glm.mod <- glm(label ~ mRNA + ORF + Fickett + Hexamer, data = subset, family = binomial(link="logit"))
#                 res <- predict(glm.mod, test.set, type = "response")
#                 res <- data.frame(Prob = res, label = test.set$label)
#
#                 res$Pred <- ifelse(res$Prob >= best.cutoff, "Coding", "NonCoding")
#                 res$label <- ifelse(res$label == 1, "Coding", "NonCoding")
#                 confusion.res <- caret::confusionMatrix(res$Pred, res$label, mode = "everything", positive = "NonCoding")
#                 performance.res <- data.frame(Sensitivity = confusion.res$byClass[1],
#                                               Specificity = confusion.res$byClass[2],
#                                               Accuracy    = confusion.res$overall[1],
#                                               F.Measure   = confusion.res$byClass[7],
#                                               Kappa       = confusion.res$overall[2])
#         })
#
#         Ave.res <- apply(performance, 1, as.numeric)
#         Ave.res <- as.data.frame(t(Ave.res))
#         Ave.res <- rowMeans(Ave.res)
#         output  <- list(Best.Cutoff = best.cutoff, Performance = Ave.res)
# }

###### Reference frequencies ######

#' Make Frequencies File for Log.Dist, Euc.Dist, and hexamer score
#'
#' @description This function is used to calculate the frequencies of lncRNAs and CDs.
#' The Frequencies file can be used to calculate Logarithm-Distance (\code{\link{compute_LogDistance}}),
#' Euclidean-Distance (\code{\link{compute_EucDistance}}), and hexamer score (\code{\link{compute_hexamerScore}}).
#'
#' NOTE: If users need to make frequencies file to build
#' new LncFinder classifier using function \code{\link{extract_features}},
#' please refer to function \code{make_frequencies}.
#'
#' @param cds.seq Coding sequences (mRNA without UTRs). Can be a FASTA file loaded
#' by \code{\link[seqinr]{seqinr-package}}.
#'
#' @param lncRNA.seq Long non-coding RNA sequences. Can be a FASTA file loaded by
#' \code{\link[seqinr]{seqinr-package}}.
#'
#' @param k An integer that indicates the sliding window size. (Default: \code{6})
#'
#' @param step Integer defaulting to \code{1} for the window step.
#'
#' @param alphabet A vector of single characters that specify the different character
#' of the sequence. (Default: \code{alphabet = c("a", "c", "g", "t")})
#'
#' @param on.orf Logical. Incomplete CDs can lead to a false shift and a
#' inaccurate hexamer frequencies. When \code{on.orf = TRUE}, the frequencies
#' will be calculated on the longest ORF. This parameter is strongly recommended to
#' set as \code{TRUE} when mRNA is used as CDs. Only available when
#' \code{alphabet = c("a", "c", "g", "t")}. (Default: \code{TRUE})
#'
#' @param ignore.illegal Logical. If \code{TRUE}, the sequences with non-nucleotide
#' characters (nucleotide characters: "a", "c", "g", "t") will be ignored when
#' calculating the frequencies. Only available when \code{alphabet = c("a", "c", "g", "t")}.
#' (Default: \code{TRUE})
#'
#' @return Returns a list which consists the frequencies of protein-coding sequences
#' and non-coding sequences.
#'
#' @author HAN Siyu
#' @details This function is used to make frequencies file for the computation of
#' Logarithm-Distance (\code{\link{compute_LogDistance}}), Euclidean-Distance
#' (\code{\link{compute_EucDistance}}),
#' and hexamer score (\code{\link{compute_hexamerScore}}).
#'
#' In order to achieve high accuracy, mRNA should not be regarded as CDs and assigned
#' to parameter \code{cds.seq}. However, CDs of some species may be insufficient
#' for calculating frequencies. In that case, mRNAs can be regarded as CDs with parameter
#' \code{on.orf = TRUE}, and the frequencies will be calculated on ORF region.
#' If \code{on.orf = TRUE}, users can set \code{step = 3} to simulate the translation process.
#'
#' @section References:
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#'
#' @seealso \code{\link{make_frequencies}},
#'          \code{\link{compute_LogDistance}},
#'          \code{\link{compute_EucDistance}},
#'          \code{\link{compute_hexamerScore}}.
#'
#' @examples
#' \dontrun{
#' Seqs <- seqinr::read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#'
#' referFreq <- make_referFreq(cds.seq = Seqs, lncRNA.seq = Seqs, k = 6, step = 1,
#'                             alphabet = c("a", "c", "g", "t"), on.orf = TRUE,
#'                             ignore.illegal = TRUE)
#' }
#'
#' @export

make_referFreq <- function(cds.seq, lncRNA.seq, k = 6, step = 1, alphabet = c("a", "c", "g", "t"),
                           on.orf = TRUE, ignore.illegal = TRUE) {

        if(on.orf & !all(alphabet %in% c("a", "t", "g", "c"))) {
                stop("Error: The cds.seq has to be DNA sequences!")
        }

        if(!all(alphabet %in% c("a", "t", "g", "c"))) {
                ignore.illegal = FALSE
        }

        message("+ Calculating frequencies. ", Sys.time(), "\n")

        if(on.orf) {
                message("- Finding ORFs")
                orf.seq <- lapply(cds.seq, function(x) {
                        orf.info <- find_ORF(x)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                cds.seq <- orf.seq[!is.na(orf.seq)]
        }

        message("- Calculating frequencies.")

        message("\n", "  Non-coding sequences")
        lnc.freq <- get.freq(lncRNA.seq, alphabet = alphabet, wordsize = k,
                             step = step, freq = F, ignore.illegal = ignore.illegal)
        message("  Completed.", "\n")

        message("  Coding sequences")
        cds.freq <- get.freq(cds.seq, alphabet = alphabet, wordsize = k,
                             step = step, freq = F, ignore.illegal = ignore.illegal)
        message("  Completed.", "\n")

        message("+ Frequencies calculation completed. ", Sys.time())

        Internal.freq <- list(ref.lnc = lnc.freq, ref.cds = cds.freq)
        Internal.freq
}

###### Logarithm-Distance ######

#' Compute Logarithm Distance
#'
#' @description This function can compute Logarithm Distance proposed by method LncFinder
#' (Han et al. 2018). Logarithm Distance can be calculated on full sequence or the longest ORF
#' region. The step and \emph{k} of the sliding window can also be customized.
#'
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
#'
#' @param label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#'
#' @param referFreq a list obtained from function \code{\link{make_referFreq}}.
#'
#' @param k An integer that indicates the sliding window size. (Default: \code{6})
#'
#' @param step Integer defaulting to \code{1} for the window step.
#'
#' @param alphabet A vector of single characters that specify the different character
#' of the sequence. (Default: \code{alphabet = c("a", "c", "g", "t")})
#'
#' @param on.ORF  Logical. If \code{TRUE}, Logarithm Distance will be calculated on
#' the longest ORF region. NOTE: If \code{TRUE}, the input has to be DNA sequences.
#' (Default: \code{FALSE})
#'
#' @param auto.full Logical. When \code{on.ORF = TRUE} but no ORF can be found,
#' if \code{auto.full = TRUE}, Logarithm Distance will be calculated on full sequences automatically;
#' if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded. Ignored when \code{on.ORF = FALSE}.
#' (Default: \code{FALSE})
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return A dataframe.
#' @author HAN Siyu
#' @details This function can compute Logarithm Distance proposed by LncFinder (HAN et al. 2018).
#' In LncFinder, two schemes are provided to calculate Logarithm Distance:
#' 1) \code{step = 3} and \code{k = 6} on the longest ORF region;
#' 2) \code{step = 1} and \code{k = 6} on full sequence.
#' Method LncFinder uses scheme 1 to extract Logarithm Distance features.
#' Using this function \code{compute_EucDistance}, both \code{step}, \code{k},
#' and calculated region (full sequence or ORF)
#' can be customized to maximize its availability.
#'
#' @section References:
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @seealso \code{\link{make_referFreq}},
#'          \code{\link{compute_EucDistance}},
#'          \code{\link{compute_hexamerScore}}.
#'
#' @examples
#' \dontrun{
#' Seqs <- seqinr::read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#'
#' referFreq <- make_referFreq(cds.seq = Seqs, lncRNA.seq = Seqs, k = 6, step = 3,
#'                             alphabet = c("a", "c", "g", "t"), on.orf = TRUE,
#'                             ignore.illegal = TRUE)
#'
#' data(demo_DNA.seq)
#' Sequences <- demo_DNA.seq
#'
#' LogDistance <- compute_LogDistance(Sequences, label = "NonCoding", referFreq = referFreq,
#'                                    k = 6, step = 3, alphabet = c("a", "c", "g", "t"),
#'                                    on.ORF = TRUE, auto.full = TRUE, parallel.cores = 2)
#' }
#'
#' @export

compute_LogDistance <- function(Sequences, label = NULL, referFreq,
                                k = 6, step = 1, alphabet = c("a", "c", "g", "t"),
                                on.ORF = FALSE, auto.full = FALSE, parallel.cores = 2) {

        if(on.ORF & !all(alphabet %in% c("a", "t", "g", "c"))) {
                stop("Error: Calculation of Logarithm Distance on ORF region can only be applied to DNA sequences!")
        }

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        message("Processing... ", Sys.time(), "\n")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = "find_ORF", envir = environment())

        if(on.ORF & !auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                Sequences <- Sequences[!is.na(Sequences)]
        }

        if(on.ORF & auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- x
                        orf
                })
        }

        message("Calculating Logarithm Distance...")
        seqFreq <- parallel::parSapply(cl, Sequences, Internal.LogDistance, k = k, step = step,
                                       alphabet = alphabet, referFreq = referFreq)
        parallel::stopCluster(cl)

        seqFreq.df <- data.frame(t(seqFreq))

        if(!is.null(label)) seqFreq.df <- cbind(label = label, seqFreq.df)
        message("\n", "Completed. ", Sys.time(), "\n")

        seqFreq.df
}

###### Euclidean-Distance ######

#' Compute Euclidean Distance
#'
#' @description This function can compute Euclidean Distance proposed by method LncFinder
#' (Han et al. 2018). Euclidean Distance can be calculated on full sequence or the longest ORF
#' region. The step and \emph{k} of the sliding window can also be customized.
#'
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
#'
#' @param label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#'
#' @param referFreq a list obtained from function \code{\link{make_referFreq}}.
#'
#' @param k An integer that indicates the sliding window size. (Default: \code{6})
#'
#' @param step Integer defaulting to \code{1} for the window step.
#'
#' @param alphabet A vector of single characters that specify the different character
#' of the sequence. (Default: \code{alphabet = c("a", "c", "g", "t")})
#'
#' @param on.ORF  Logical. If \code{TRUE}, Euclidean Distance will be calculated on
#' the longest ORF region. NOTE: If \code{TRUE}, the input has to be DNA sequences.
#' (Default: \code{FALSE})
#'
#' @param auto.full Logical. When \code{on.ORF = TRUE} but no ORF can be found,
#' if \code{auto.full = TRUE}, Euclidean Distance will be calculated on full sequences automatically;
#' if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded. Ignored when \code{on.ORF = FALSE}.
#' (Default: \code{FALSE})
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return A dataframe.
#' @author HAN Siyu
#' @details This function can compute Euclidean Distance proposed by LncFinder (HAN et al. 2018).
#' In LncFinder, two schemes are provided to calculate Euclidean Distance:
#' 1) \code{step = 3} and \code{k = 6} on the longest ORF region;
#' 2) \code{step = 1} and \code{k = 6} on full sequence.
#' Using this function \code{compute_EucDistance}, both \code{step}, \code{k},
#' and calculated region (full sequence or ORF)
#' can be customized to maximize its availability.
#'
#' @section References:
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @seealso \code{\link{make_referFreq}},
#'          \code{\link{compute_LogDistance}},
#'          \code{\link{compute_hexamerScore}}.
#'
#' @examples
#' \dontrun{
#' Seqs <- seqinr::read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#'
#' referFreq <- make_referFreq(cds.seq = Seqs, lncRNA.seq = Seqs, k = 6, step = 3,
#'                             alphabet = c("a", "c", "g", "t"), on.orf = TRUE,
#'                             ignore.illegal = TRUE)
#'
#' data(demo_DNA.seq)
#' Sequences <- demo_DNA.seq
#'
#' EucDistance <- compute_EucDistance(Sequences, label = "NonCoding", referFreq = referFreq,
#'                                    k = 6, step = 3, alphabet = c("a", "c", "g", "t"),
#'                                    on.ORF = TRUE, auto.full = TRUE, parallel.cores = 2)
#' }
#'
#' @export

compute_EucDistance <- function(Sequences, label = NULL, referFreq,
                                k = 6, step = 1, alphabet = c("a", "c", "g", "t"),
                                on.ORF = FALSE, auto.full = FALSE, parallel.cores = 2) {

        if(on.ORF & !all(alphabet %in% c("a", "t", "g", "c"))) {
                stop("Error: Calculation of Euclidean Distance on ORF region can only be applied to DNA sequences!")
        }

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        message("Processing... ", Sys.time(), "\n")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = "find_ORF", envir = environment())

        if(on.ORF & !auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                Sequences <- Sequences[!is.na(Sequences)]
        }

        if(on.ORF & auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- x
                        orf
                })
        }

        message("Calculating Euclidean Distance...")
        seqFreq <- parallel::parSapply(cl, Sequences, Internal.EucDistance, k = k, step = step,
                                       alphabet = alphabet, referFreq = referFreq)
        parallel::stopCluster(cl)

        seqFreq.df <- data.frame(t(seqFreq))

        if(!is.null(label)) seqFreq.df <- cbind(label = label, seqFreq.df)
        message("\n", "Completed. ", Sys.time(), "\n")

        seqFreq.df
}

###### Hexamer Score ######

#' Compute Hexamer Score
#'
#' @description This function can compute hexamer score proposed by method CPAT
#' (Wang et al. 2013). Hexamer score can be calculated on full sequence or the longest ORF
#' region. The step and \emph{k} of the sliding window can also be customized.
#'
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
#'
#' @param label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#'
#' @param referFreq A list obtained from function \code{\link{make_referFreq}}.
#'
#' @param k An integer that indicates the sliding window size. (Default: \code{6})
#'
#' @param step Integer defaulting to \code{1} for the window step.
#'
#' @param alphabet A vector of single characters that specify the different character
#' of the sequence. (Default: \code{alphabet = c("a", "c", "g", "t")})
#'
#' @param on.ORF Logical. If \code{TRUE}, hexamer score will be calculated on
#' the longest ORF region. NOTE: If \code{TRUE}, the input has to be DNA sequences.
#' (Default: \code{FALSE})
#'
#' @param auto.full Logical. When \code{on.ORF = TRUE} but no ORF can be found,
#' if \code{auto.full = TRUE},  hexamer score will be calculated on full sequences automatically;
#' if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded. Ignored when \code{on.ORF = FALSE}.
#' (Default: \code{FALSE})
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return A dataframe.
#' @author HAN Siyu
#' @details This function can compute hexamer score proposed by CPAT (Wang et al. 2013).
#' In CPAT, hexamer score is calculated on the longest ORF region, and the step of the
#' sliding window is 3 (i.e. \code{step = 3}). Hexamer means six adjoining bases, thus
#' \code{k = 6}. But in function \code{compute_hexamerScore}, both \code{step}, \code{k},
#' and calculated region (full sequence or ORF)
#' can be customized to maximize its availability.
#'
#' @section References:
#' Liguo Wang, Hyun Jung Park, Surendra Dasari, Shengqin Wang, JeanPierre Kocher, & Wei Li.
#' CPAT: coding-potential assessment tool using an alignment-free logistic regression model.
#' \emph{Nucleic Acids Research}, 2013, 41(6):e74-e74.
#'
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @seealso \code{\link{make_referFreq}},
#'          \code{\link{compute_LogDistance}},
#'          \code{\link{compute_EucDistance}}.
#'
#' @examples
#' \dontrun{
#' Seqs <- seqinr::read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#'
#' referFreq <- make_referFreq(cds.seq = Seqs, lncRNA.seq = Seqs, k = 6, step = 1,
#'                             alphabet = c("a", "c", "g", "t"), on.orf = TRUE,
#'                             ignore.illegal = TRUE)
#'
#' data(demo_DNA.seq)
#' Sequences <- demo_DNA.seq
#'
#' hexamerScore <- compute_hexamerScore(Sequences, label = "NonCoding", referFreq = referFreq,
#'                                      k = 6, step = 1, alphabet = c("a", "c", "g", "t"),
#'                                      on.ORF = TRUE, auto.full = TRUE, parallel.cores = 2)
#' }
#'
#' @export

compute_hexamerScore <- function(Sequences, label = NULL, referFreq,
                                 k = 6, step = 1, alphabet = c("a", "c", "g", "t"),
                                 on.ORF = FALSE, auto.full = FALSE, parallel.cores = 2) {

        if(on.ORF & !all(alphabet %in% c("a", "t", "g", "c"))) {
                stop("Error: Calculation of hexamer score on ORF region can only be applied to DNA sequences!")
        }

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        message("Processing... ", Sys.time(), "\n")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = "find_ORF", envir = environment())

        if(on.ORF & !auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                Sequences <- Sequences[!is.na(Sequences)]
        }

        if(on.ORF & auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- x
                        orf
                })
        }

        message("Calculating hexamer score...")
        seqFreq <- parallel::parSapply(cl, Sequences, Internal.hexamerScore, k = k, step = step,
                                       alphabet = alphabet, referFreq = referFreq)
        parallel::stopCluster(cl)

        seqFreq.df <- data.frame(Hexamer.Score = seqFreq)

        if(!is.null(label)) seqFreq.df <- cbind(label = label, seqFreq.df)
        message("\n", "Completed. ", Sys.time(), "\n")

        seqFreq.df
}

###### k-mer scheme ######

#' Compute \emph{k}-mer Features
#'
#' @description This function can calculate the \emph{k}-mer frequencies of the sequences.
#'
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
#'
#' @param label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#'
#' @param k An integer that indicates the sliding window size. (Default: \code{1:5})
#'
#' @param step Integer defaulting to \code{1} for the window step.
#'
#' @param freq Logical. If TRUE, the frequencies of different patterns are returned
#' instead of counts. (Default: \code{TRUE})
#'
#' @param improved.mode Logical. If TRUE, the frequencies will be normalized using
#' the method proposed by PLEK (Li et al. 2014).
#' Ignored if \code{freq = FALSE}. (Default: \code{FALSE})
#'
#' @param alphabet A vector of single characters that specify the different character
#' of the sequence. (Default: \code{alphabet = c("a", "c", "g", "t")})
#'
#' @param on.ORF Logical. If \code{TRUE}, the \emph{k}-mer frequencies will be calculated on
#' the longest ORF region. NOTE: If \code{TRUE}, the sequences have to be DNA.
#' (Default: \code{FALSE})
#'
#' @param auto.full Logical. When \code{on.ORF = TRUE} but no ORF can be found,
#' if \code{auto.full = TRUE}, the \emph{k}-mer
#' frequencies will be calculated on the full sequence automatically;
#' if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded.
#' Ignored when \code{on.ORF = FALSE}. (Default: \code{FALSE})
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return A dataframe.
#' @author HAN Siyu
#' @details
#' This function can extract \emph{k}-mer features. \code{k} and \code{step} can be customized.
#' The count (\code{freq = FALSE}) or frequencies (\code{freq = TRUE}) of different patterns can be returned.
#' If \code{freq = TRUE}, \code{improved.mode} is available. The improved mode is proposed by method PLEK.
#' (Ref: Li et al. 2014)
#'
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' kmer_res1 <- compute_kmer(Seqs, k = 1:5, step = 1, freq = TRUE, improved.mode = FALSE)
#'
#' kmer_res2 <- compute_kmer(Seqs, k = 1:5, step = 3, freq = TRUE,
#'                           improved.mode = TRUE, on.ORF = TRUE, auto.full = TRUE)
#' }
#'
#' @export

compute_kmer <- function(Sequences, label = NULL, k = 1:5, step = 1, freq = TRUE,
                         improved.mode = FALSE, alphabet = c("a", "c", "g", "t"),
                         on.ORF = FALSE, auto.full = FALSE, parallel.cores = 2) {

        if(on.ORF & !all(alphabet %in% c("a", "t", "g", "c"))) {
                stop("Error: Calculation of k-mer frequencies on ORF region can only be applied to DNA sequences!")
        }

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        message("Processing... ", Sys.time(), "\n")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = "find_ORF", envir = environment())

        if(on.ORF & !auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                Sequences <- Sequences[!is.na(Sequences)]
        }

        if(on.ORF & auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- x
                        orf
                })
        }

        for (i in k) {
                message("Calculating ", i, "-mer frequencies...")
                i.freq <- parallel::parSapply(cl, Sequences, seqinr::count, wordsize = i,
                                              by = step, freq = freq, alphabet = alphabet)
                if(improved.mode & freq) {
                        weight <- 1 / (length(alphabet) ^ (max(k) - i))
                        i.freq <- i.freq * weight
                }
                i.freq <- data.frame(t(i.freq))
                if(i == k[1]) {
                        k.freq <- i.freq
                } else {
                        k.freq <- cbind(k.freq, i.freq)
                }
        }
        parallel::stopCluster(cl)

        if(!is.null(label)) k.freq <- cbind(label = label, k.freq)
        message("\n", "Completed. ", Sys.time(), "\n")
        k.freq
}

###### Fickett Score ######

#' Compute Fickett TESTCODE Score
#'
#' @description This function can compute Fickett TESTCODE score of DNA sequences proposed by James W.Fickett
#' (Fickett JW. 1982). Fickett TESTCODE score can be calculated on full sequence or the longest ORF
#' region.
#'
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
#'
#' @param label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#'
#' @param on.ORF Logical. If \code{TRUE}, Fickett TESTCODE score will be calculated on
#' the longest ORF region.
#'
#' @param auto.full Logical. When \code{on.ORF = TRUE} but no ORF can be found,
#' if \code{auto.full = TRUE},  Fickett TESTCODE score will be calculated on full sequences automatically;
#' if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded.
#' Ignored when \code{on.ORF = FALSE}. (Default: \code{FALSE})
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return A dataframe.
#' @author HAN Siyu
#' @details This function can compute Fickett TESTCODE score proposed by James W.Fickett (Fickett JW. 1982).
#' Fickett TESTCODE score is selected as feature by method CPAT (Wang et al. 2013) and CPC2 (Kang et al. 2017).
#' In CPAT, Fickett TESTCODE score is calculated on the longest ORF region, but CPC2 calculates the score
#' on full sequence. This function \code{compute_FickettScore} improves the CPAT's code
#' and is capable of computing the score on the longest ORF region as well as full sequence.
#'
#' @section References:
#' James W.Fickett.
#' Recognition of protein coding regions in DNA sequences.
#' \emph{Nucleic Acids Research}, 1982, 10(17):5303-5318.
#'
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
#'
#' Liguo Wang, Hyun Jung Park, Surendra Dasari, Shengqin Wang, JeanPierre Kocher & Wei Li.
#' CPAT: coding-potential assessment tool using an alignment-free logistic regression model.
#' \emph{Nucleic Acids Research}, 2013, 41(6):e74-e74.
#'
#' Yu-Jian Kang, De-Chang Yang, Lei Kong, Mei Hou, Yu-Qi Meng, Liping Wei & Ge Gao.
#' CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features.
#' \emph{Nucleic Acids Research}, 2017, 45(W1):W12-W16.
#'
#' @importFrom seqinr count
#' @importFrom seqinr s2c
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' FickettScore <- compute_FickettScore(Seqs, label = NULL, on.ORF = TRUE,
#'                                      auto.full = TRUE, parallel.cores = 2)
#' }
#' @export

compute_FickettScore <- function(Sequences, label = NULL, on.ORF = FALSE,
                                 auto.full = FALSE, parallel.cores = 2) {

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        message("Processing... ", Sys.time(), "\n")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = c("find_ORF", "Internal.convertProb"), envir = environment())

        if(on.ORF & !auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                Sequences <- Sequences[!is.na(Sequences)]
        }

        if(on.ORF & auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- x
                        orf
                })
        }

        positionProb = list(a = c(0.94,0.68,0.84,0.93,0.58,0.68,0.45,0.34,0.20,0.22),
                            c = c(0.80,0.70,0.70,0.81,0.66,0.48,0.51,0.33,0.30,0.23),
                            g = c(0.90,0.88,0.74,0.64,0.53,0.48,0.27,0.16,0.08,0.08),
                            t = c(0.97,0.97,0.91,0.68,0.69,0.44,0.54,0.20,0.09,0.09))
        positionWeight = c(a = 0.26, c = 0.18, g = 0.31, t = 0.33)
        positionParam = c(1.9, 1.8, 1.7, 1.6, 1.5 ,1.4, 1.3, 1.2, 1.1, 0)

        contentProb = list(a = c(0.28,0.49,0.44,0.55,0.62,0.49,0.67,0.65,0.81,0.21),
                           c = c(0.82,0.64,0.51,0.64,0.59,0.59,0.43,0.44,0.39,0.31),
                           g = c(0.40,0.54,0.47,0.64,0.64,0.73,0.41,0.41,0.33,0.29),
                           t = c(0.28,0.24,0.39,0.40,0.55,0.75,0.56,0.69,0.51,0.58))
        contentWeight = c(a = 0.11, c = 0.12, g = 0.15, t = 0.14)
        contentParam = c(0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0)

        message("Calculating Fickett TESTCODE score...")
        FickettScore.res <- parallel::parSapply(cl, Sequences, Internal.FickettScore,
                                                .positionProb = positionProb, .positionWeight = positionWeight,
                                                .positionParam = positionParam, .contentProb = contentProb,
                                                .contentWeight = contentWeight, .contentParam = contentParam)
        parallel::stopCluster(cl)

        FickettScore.df <- data.frame(Fickett.Score = FickettScore.res)

        if(!is.null(label)) FickettScore.df <- cbind(label = label, FickettScore.df)
        message("\n", "Completed. ", Sys.time(), "\n")
        FickettScore.df
}

###### GC Content ######

#' Calculate GC content
#' @description This function can GC content of the input sequences.
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
#'
#' @param label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#'
#' @param on.ORF Logical. If \code{TRUE}, GC content will be calculated on
#' the longest ORF region.
#'
#' @param auto.full Logical. When \code{on.ORF = TRUE} but no ORF can be found,
#' if \code{auto.full = TRUE},  GC content will be calculated on full sequences automatically;
#' if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded.
#' Ignored when \code{on.ORF = FALSE}. (Default: \code{FALSE})
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return A dataframe.
#'
#' @author HAN Siyu
#'
#' @details This function can basically compute GC content of DNA sequences:
#' GC content = (nc + ng) / (na + nc + ng + nt).
#' The function will ignored the ambiguous bases.
#'
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @seealso \code{\link[seqinr]{GC}} (package "\code{\link[seqinr]{seqinr-package}}")
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' gcContent <- compute_GC(Seqs, label = "NonCoding",on.ORF = TRUE,
#'                         auto.full = TRUE, parallel.cores = 2)
#' }
#' @export
#'

compute_GC <- function(Sequences, label = NULL, on.ORF = FALSE,
                       auto.full = FALSE, parallel.cores = 2) {

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        message("Processing... ", Sys.time(), "\n")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = c("find_ORF"), envir = environment())

        if(on.ORF & !auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                Sequences <- Sequences[!is.na(Sequences)]
        }

        if(on.ORF & auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- x
                        orf
                })
        }

        message("Calculating GC content...")
        GC.content <- parallel::parSapply(cl, Sequences, function(x) {
                num.a <- sum(x == "a")
                num.c <- sum(x == "c")
                num.g <- sum(x == "g")
                num.t <- sum(x == "t")

                gc.content <- (num.g + num.c) / (num.a + num.c + num.g + num.t)
        })
        parallel::stopCluster(cl)

        GC.content  <- as.data.frame(GC.content, stringsAsFactors = FALSE)

        if(!is.null(label)) GC.content <- cbind(label = label, GC.content)

        message("\n", "Completed. ", Sys.time(), "\n")

        GC.content
}

###### EIIP Features ######

#' Extract the EIIP-derived features
#'
#' @description This function can extract EIIP-derived features proposed by Han et al (2018).
#'
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
#' @param label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#' @param spectrum.percent Numeric specifying the percentage of the sorted power spectrum that be
#' used to calculate the quantile-based features. For example, if \code{spectrum.percent = 0.1},
#' the top 10\% percent of the sorted power spectrum will be used to compute the quantiles.
#' @param quantile.probs Numeric. The probabilities with values in [0,1].
#'
#' @return A dataframe including the EIIP-derived features.
#'
#' @author HAN Siyu
#'
#' @details The function \code{compute_EIIP} can extract EIIP (electron-ion interaction pseudo-potential) features including:
#' signal at 1/3 position (\code{Signal.Peak}), average power (\code{Average.Power}), signal to noise ratio (\code{SNR}),
#' and quantile-based features of one specified percentage of the sorted power spectrum
#' (e.g. \code{0\%}, \code{20\%}, \code{40\%}, \code{60\%}, \code{70\%},
#' \code{100\%} when \code{quantile.probs = seq(0, 1, 0.2)} and \code{spectrum.percent =} \code{0.1}).
#'
#' In method LncFinder, EIIP features includes \code{Signal.Peak}, \code{SNR}, 0\% (\code{Signal.Min}),
#' 25\% (\code{Singal.Q1}, 50\% \code{Signal.Q2}), and 75\% (\code{Signal.Max}) of the top 10\% sorted
#' power spectrum, i.e. \code{quantile.prob} \code{= seq(0, 1, 0.25)} and \code{spectrum.percent = 0.1}.
#'
#' @importFrom stats fft
#' @importFrom stats quantile
#'
#' @seealso \code{\link{extract_features}}
#'
#' @section References:
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
#'
#' Achuthsankar S Nair & Sivarama Pillai Sreenadhan.
#' A coding measure scheme employing electron-ion interaction pseudopotential (EIIP).
#' \emph{Bioinformation}, 2006, 1(6):197-202.
#'
#' @examples
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' EIIP_res <- compute_EIIP(Seqs, label = "NonCoding", spectrum.percent = 0.25,
#'                          quantile.probs = seq(0, 1, 0.25))
#'
#' @export

compute_EIIP <- function(Sequences, label = NULL, spectrum.percent = 0.1,
                         quantile.probs = seq(0, 1, 0.25)) {
        seqs.EIIP <- sapply(Sequences, get_EIIP, percent.length = spectrum.percent, probs = quantile.probs)

        seqs.EIIP  <- as.data.frame(t(seqs.EIIP), stringsAsFactors = FALSE)

        if(!is.null(label)) seqs.EIIP <- cbind(label = label, seqs.EIIP)

        seqs.EIIP
}

###### pI Features ######

#' Compute Theoretical Isoelectric Point
#'
#' @description This function is basically a wrapper for function \code{\link[seqinr]{computePI}}.
#' This function translate DNA sequence into protein, and compute the theoretical isoelectric point
#' (pI) of this protein.
#'
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
#'
#' @param label Optional. String. Indicate the label of the sequences such as
#' "NonCoding", "Coding".
#'
#' @param on.ORF Logical. If \code{TRUE}, pI will be calculated on
#' the longest ORF region. NOTE: If \code{TRUE}, the input has to be DNA sequences.
#' (Default: \code{FALSE})
#'
#' @param auto.full Logical. When \code{on.ORF = TRUE} but no ORF can be found,
#' if \code{auto.full = TRUE},  pI will be calculated on full sequences automatically;
#' if \code{auto.full} is \code{FALSE}, the sequences that have no ORF will be discarded.
#' Ignored when \code{on.ORF = FALSE}. (Default: \code{FALSE})
#'
#' @param ambiguous.base If \code{TRUE}, ambiguous bases are taken into account when
#' translating DNA sequences into proteins.
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return A dataframe.
#' @author HAN Siyu
#' @details This function can compute the pI of DNA sequences. Method CPC2 (Kang et al. 2017) uses this
#' feature to identify lncRNAs, and this feature is evaluated in the article LncFinder (Han et al. 2018).
#'
#' Using this function, the theoretical pI can be computed on full sequence or the longest ORF region.
#' In CPC2, pI is calculated on ORF region.
#'
#' @importFrom seqinr translate
#' @importFrom seqinr computePI
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Sequences <- demo_DNA.seq
#'
#' pI_res <- compute_pI(Sequences, on.ORF = TRUE, auto.full = FALSE, ambiguous.base = FALSE)
#' }
#' @export

compute_pI <- function(Sequences, label = NULL, on.ORF = FALSE, auto.full = FALSE,
                       ambiguous.base = FALSE, parallel.cores = 2) {

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        message("Processing... ", Sys.time(), "\n")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = "find_ORF", envir = environment())

        if(on.ORF & !auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                Sequences <- Sequences[!is.na(Sequences)]
        }

        if(on.ORF & auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- x
                        orf
                })
        }

        message("Calculating theoretical isoelectric point...")
        AA.seq  <- parallel::parLapply(cl, Sequences, seqinr::translate, ambiguous = ambiguous.base)
        pI.val  <- parallel::parSapply(cl, AA.seq,  seqinr::computePI)

        parallel::stopCluster(cl)

        pI      <- data.frame(pI = pI.val, stringsAsFactors = FALSE)

        if(!is.null(label)) pI <- cbind(label = label, pI)
        message("\n", "Completed. ", Sys.time(), "\n")

        pI
}

###### SVM k-fold CV ######

#' \emph{k}-fold Cross Validation for SVM
#' @description This function conduct \emph{k}-fold Cross Validation for SVM.
#'
#' @param dataset The dataset obtained from function \code{\link{extract_features}}.
#' Or datasets used to build the classifier.
#'
#' @param label.col integer specifying the column number of the label. (Default: \code{1})
#'
#' @param positive.class Character. Indicate the positive class of the dataset.
#' (Default: \code{NonCoding}) The value of this parameter should be identical to
#' one of the classes of the response vectors.
#'
#' @param folds.num Integer. Specify the number of folds for cross-validation.
#' (Default: \code{10})
#'
#' @param seed Integer. Used to set the seed for cross-validation. (Default: \code{1})
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}, users can set as \code{-1} to run
#' this function with all cores. If the number of \code{parallel.cores} is more
#' than the \code{folds.num} (number of the folds for cross-validation), the
#' number of \code{parallel.cores} will be set as \code{folds.num} automatically.
#'
#' @param ... additional parameters for function \code{\link[e1071]{svm}}.
#'
#' @return Returns the optimal parameters when \code{return.model = FALSE}.
#' Or returns the best model when \code{return.model = TRUE}.
#'
#' @author HAN Siyu
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
#' @seealso \code{\link{extract_features}}, \code{\link{svm_tune}}.
#' @examples
#' \dontrun{
#' data(demo_dataset)
#' my_dataset <- demo_dataset
#'
#' cv_res <- svm_cv(my_dataset, folds.num = 4, seed = 1,
#'                  parallel.core = 2, cost = 3, kernel = "radial", gamma = 0.5)
#'
#' ### Users can set return.model = TRUE to return the best model.
#' }
#' @export

svm_cv <- function(dataset, label.col = 1, positive.class = NULL,
                   folds.num = 10, seed = 1, parallel.cores = 2, ...) {

        names(dataset)[[label.col]] <- "label"
        set.seed(seed)
        folds <- caret::createFolds(dataset$label, k = folds.num, returnTrain = TRUE)

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        parallel.cores <- ifelse(parallel.cores > folds.num, folds.num, parallel.cores)

        if(parallel.cores == 2 & folds.num > parallel.cores) {
                message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        }

        cl <- parallel::makeCluster(parallel.cores, outfile = "")

        if(is.null(positive.class)) {
                positive.class <- as.character(dataset$label[[1]])
                parallel::clusterExport(cl, varlist = c("dataset", "positive.class"), envir = environment())

        } else {
                parallel::clusterExport(cl, varlist = c("dataset"), envir = environment())
        }
        list(...)

        perf.res <- parallel::parSapply(cl, folds, fold.res, dataset = dataset,
                                        positive.class = positive.class, ...)
        parallel::stopCluster(cl)
        Ave.res2     <- apply(perf.res, 1, as.numeric)
        Ave.res2     <- as.data.frame(t(Ave.res2))
        Ave.res2$Ave.Res <- rowMeans(Ave.res2)
        message("Average Result:")
        print(Ave.res2[ncol(Ave.res2)])
        Ave.res2
        # perf.res
}

###### svm_tune ######

#' Parameter Tuning of SVM
#' @description This function conduct the parameter tuning of SVM. Parameters
#' gamma and cost can be tuned using grid search.
#'
#' @param dataset The dataset obtained from function \code{\link{extract_features}}.
#' Or datasets used to build the classifier.
#'
#' @param label.col integer specifying the column number of the label. (Default: \code{1})
#'
#' @param positive.class Character. Indicate the positive class of the dataset.
#' (Default: \code{NonCoding}) The value of this parameter should be identical to
#' one of the classes of the response vectors.
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
#' @return Returns the optimal parameters when \code{return.model = FALSE}.#'
#' Or returns the best model when \code{return.model = TRUE}.
#'
#' @author HAN Siyu
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
#' @seealso \code{\link{extract_features}}, \code{\link{svm_cv}}.
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' positive_data <- extract_features(Seqs[1:5], label = "NonCoding",
#'                                   SS.features = FALSE, format = "DNA",
#'                                   frequencies.file = "human",
#'                                   parallel.cores = 2)
#'
#' negative_data <- extract_features(Seqs[6:10], label = "Coding",
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
#' optimal_parameter <- svm_tune(my_dataset, positive.class = "NonCoding",
#'                               folds.num = 2, seed = 1,
#'                               gamma.range = (2 ^ seq(-5, 0, 2)),
#'                               cost.range = c(1, 8, 16),
#'                               return.model = FALSE, parallel.core = 2)
#'
#' ### Users can set return.model = TRUE to return the best model.
#' }
#' @export

svm_tune <- function(dataset, label.col = 1, positive.class = "NonCoding",
                     folds.num = 10, seed = 1,
                     gamma.range = (2 ^ seq(-5, 0, 1)),
                     cost.range = c(1, 4, 8, 16, 24, 32),
                     return.model = TRUE, parallel.cores = 2){

        names(dataset)[[label.col]] <- "label"

        set.seed(seed)
        folds <- caret::createFolds(dataset$label, k = folds.num, returnTrain = TRUE)

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        parallel.cores <- ifelse(parallel.cores > folds.num, folds.num, parallel.cores)

        if(parallel.cores == 2 & folds.num > parallel.cores) {
                message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        }

        cl <- parallel::makeCluster(parallel.cores)

        res <- c()
        message("+ SVM.tune processing.")
        for(g in gamma.range){
                for(C in cost.range){
                        parallel::clusterExport(cl, varlist = c("g", "C", "dataset", "positive.class"),
                                                envir = environment())
                        message("- gamma = ", g, ", Cost = ", C)
                        perf.res    <- parallel::parSapply(cl, folds, function(x, gamma = g, cost = C,
                                                                               positive.label = positive.class) {
                                subset <- dataset[x, ]
                                test.set  <- dataset[-x, ]
                                svm.mod <- e1071::svm(label ~ ., data = subset, scale = TRUE, probability = TRUE,
                                                      kernel = "radial", gamma = gamma, cost = cost)

                                res <- stats::predict(svm.mod, test.set, probability = TRUE)

                                confusion.res <- caret::confusionMatrix(data.frame(res)$res, test.set$label,
                                                                        positive = positive.label, mode = "everything")
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
                message("\n", "+ Training the model on the Full Dataset.")
                svm.mod  <- e1071::svm(label ~ ., data = dataset, scale = TRUE,
                                       probability = TRUE, kernel = "radial",
                                       gamma = best.gamma,
                                       cost  = best.cost)
                message("\n", "+ SVM.tune completed.")
                return(svm.mod)
        } else {
                message("\n", "+ SVM.tune completed.")
                return(list(Best.Parameters = best.parameters, Result = res))
        }
}

###### find_orfs ######

#' Find ORFs
#' @description This function can find all the ORFs in one sequence.
#'
#' @param OneSeq Is one sequence. Can be a FASTA file read by package "seqinr"
#' \code{\link[seqinr]{seqinr-package}} or just a string.
#'
#' @param reverse.strand Logical. Whether find ORF on the reverse strand. Default: \code{FALSE}
#'
#' @param max.only Logical. If \code{TRUE}, only the longest ORF will be returned.  Default: \code{TRUE}
#'
#' @return If \code{max.only = TRUE}, the function returns a list which consists the ORF region (\code{ORF.Max.Seq}),
#' length (\code{ORF.Max.Len}) and coverage (\code{ORF.Max.Cov}) of the longest ORF.
#' If \code{max.only = FALSE}, the function returns a dataframe which consists all the ORF sequences.
#' @author HAN Siyu
#' @details This function can extract ORFs of one sequence. It returns
#' ORF region, length and coverage of the longest ORF when \code{max.only = TRUE} or
#' ORF region, start position, end position, length and coverage of all the ORFs when
#' \code{max.only = FALSE}. Coverage is the the ratio
#' of the ORF to transcript length. If \code{reverse.strand = TRUE}, ORF will also be
#' found on reverse strand.
#' @importFrom seqinr getSequence
#' @importFrom seqinr comp
#' @importFrom seqinr s2c
#' @examples
#' ### For one sequence:
#' OneSeq <- c("cccatgcccagctagtaagcttagcc")
#' orf.info_1 <- find_orfs(OneSeq, reverse.strand = TRUE, max.only = FALSE)
#'
#' ### For a FASTA file contains several sequences:
#' \dontrun{
#' ### Use "read.fasta" function of package "seqinr" to read a FASTA file:
#' Seqs <- seqinr::read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#' }
#'
#' ### Or just try to use our data "demo_DNA.seq"
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' ### Use apply function to find the longest ORF:
#' orf.info_2 <- sapply(Seqs, find_orfs, reverse.strand = FALSE, max.only = FALSE)
#' @export

find_orfs <- function(OneSeq, reverse.strand = FALSE, max.only = TRUE) {

        orf.info <- find_ORF(OneSeq, max.only = max.only)
        if(reverse.strand) {
                OneSeq <- unlist(seqinr::getSequence(OneSeq, as.string = TRUE))
                OneSeq <- gsub("\n", "", OneSeq)
                OneSeq <- seqinr::s2c(OneSeq)
                Seq.reverse <- rev(seqinr::comp(OneSeq))
                orf.reverse <- find_ORF(Seq.reverse, max.only = max.only)
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
#' @param Sequences A FASTA file loaded by function  \code{\link[seqinr]{read.fasta}} of
#' \code{\link[seqinr]{seqinr-package}}.
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
#' @author HAN Siyu
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
#' If users have their own SS data, users can use function \code{\link{read_SS}} to load
#' them, instead of obtaining from RNAfold.
#'
#' @importFrom seqinr getSequence
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @seealso \code{\link{read_SS}}
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
#' SS.seq_1 <- run_RNAfold(Seqs[1:2], RNAfold.path = RNAfold.path, parallel.cores = 2)
#'
#' ### For UNIX/Linux, "RNAfold.path" can be just defined as "RNAfold" as default:
#' SS.seq_2 <- run_RNAfold(Seqs, RNAfold.path = "RNAfold", parallel.cores = 2)
#' }
#' @export

run_RNAfold <- function(Sequences, RNAfold.path = "RNAfold",
                        parallel.cores = 2){

        Seqs.validate <- Sequences[which(lengths(Sequences) < 30000)]

        if(length(Seqs.validate) < length(Sequences)){
                message("Due to the limitation of RNAfold,")
                message("Sequences with length more than 30000 nt will be omitted.")
                message(length(Sequences) - length(Seqs.validate), " sequences have been removed.", "\n")
                Sequences <- Seqs.validate
        }

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores, outfile = "")

        message("\n", "Sequences Number: ", length(Sequences), "\n")
        message("Processing...", "\n")
        index <- 1
        info <- paste(ceiling(length(Sequences) / parallel.cores), ",", sep = "")
        parallel::clusterExport(cl, varlist = c("info", "index"), envir = environment())
        sec.seq <- parallel::parSapply(cl, Sequences, secondary_seq, info = info,
                                       RNAfold.path = RNAfold.path)
        parallel::stopCluster(cl)
        sec.seq <- data.frame(sec.seq, stringsAsFactors = FALSE)
        message("\n", "Completed.", "\n")
        sec.seq
}

###### read_SS ######

#' Read Secondary Structure Information
#' @description This function can read secondary structure information from your
#' own file instead of obtaining from function \code{\link{run_RNAfold}}. This function
#' will be useful if users have had secondary structure sequences (Dot-Bracket Notation).
#'
#' @param oneFile.loc String. The location of your sequence file. This file should contains
#' one (and only one) RNA sequence and its secondary structure sequence in Dot-Bracket Notation.
#' This parameter needs to be defined only when \code{separateFile = FALSE}. See Details for more
#' information.
#'
#' @param seqRNA.loc String. The location of your RNA sequences file (FASTA format). If your
#' RNA sequences and secondary structure sequences are in two files, you need to define the
#' locations of two files respectively. And the files with multiple sequences are supported
#' for this option. This parameter needs to be defined only when \code{separateFile} is \code{TRUE}.
#' Location of secondary structure sequences file is also needed (parameter \code{seqSS.loc}).
#' See Details for more information.
#'
#' @param seqSS.loc String. The location of your secondary structure sequences file (FASTA format).
#'
#' @param separateFile Logical. Your RNA sequence(s) and secondary structure sequence(s) are in
#' separate files? If \code{separateFile = FALSE}, your file should have one (and only one) RNA
#' sequence and its secondary structure sequence. No limit when \code{separateFile = TRUE}.
#'
#' @param withMFE Logical. Whether MFE is provided at the end of secondary structure sequence.
#' If \code{withMFE = TRUE}, MFE will be extracted. The format should be in accordance with
#' the output format of RNAfold.
#'
#' @return A dataframe. The first row is RNA sequence, the second row is Dot-Bracket Notation of
#' secondary structure sequence, the third row is MFE (if MFE is provided).
#'
#' @author HAN Siyu
#'
#' @details When users want to predict sequences with secondary structure features, users may have
#' had their own secondary structure sequences. With this function, users can read SS information
#' from their files. Two kind of files are supported: RNA sequence and SS sequence in one file
#' \code{separateFile} is \code{FALSE} or in separate files \code{separateFile = TRUE}.
#'
#' \code{separateFile = FALSE} is used for secondary structure that obtained from some popular
#' programs, such as RNAfold. In this case, the output file only contains one RNA sequence and
#' its SS. Besides, this file only have two rows: RNA sequence and its SS sequences. Thus, this
#' option is more favorable when the file only have one sequence and the sequence are in accordance
#' with the output format of RNAfold.
#'
#' If users obtained the SS sequence from experiments, RNA sequence and SS sequence may be in two
#' files. In this case, users can select \code{separateFile = TRUE}. Two files should be in FASTA
#' format and one file can have multiple sequences. The sequences in two files should have the same
#' order. If your data are obtained from experiments or other sources, it is highly recommended
#' that users should build new model with this data, since the SS sequences of pre-built model are
#' obtained for RNAfold and may have many differences with experimental data.
#'
#' @importFrom seqinr read.fasta
#' @importFrom seqinr getSequence
#' @seealso \code{\link{run_RNAfold}}
#' @examples
#' \dontrun{
#' ### Load sequence data
#' data("demo_DNA.seq")
#' Seqs <- demo_DNA.seq[1:4]
#' ### Convert sequences from vector to string.
#' Seqs <- sapply(Seqs, seqinr::getSequence, as.string = TRUE)
#' ### Write a fasta file.
#' seqinr::write.fasta(Seqs, names = names(Seqs), file.out = "tmp.RNA.fa", as.string = TRUE)
#'
#' ### For Windows system: (Your path of RNAfold.)
#' RNAfold.path <- '"E:/Program Files/ViennaRNA/RNAfold.exe"'
#' ### Define the parameters of RNAfold. See documents of RNAfold for more information.
#' RNAfold.command <- paste(RNAfold.path, " --noPS -i tmp.RNA.fa -o output", sep = "")
#' ### Run RNAfold and output four result files.
#' system(RNAfold.command)
#'
#' ### Read secondary structure information for one file.
#' result_1 <- read_SS(oneFile.loc = "output_ENST00000510062.1.fold",
#'                     separateFile = FALSE, withMFE = TRUE)
#' ### Read secondary sturcture sequences for multiple files.
#' filePath <- dir(pattern = ".fold")
#' result_2 <- sapply(filePath, read_SS, separateFile = FALSE, withMFE = TRUE)
#' result_2 <- as.data.frame(result_2)
#' }
#'
#' @export
#'

read_SS <- function(oneFile.loc, seqRNA.loc, seqSS.loc, separateFile = TRUE, withMFE = TRUE) {
        if(separateFile) {
                RNA.seq <- seqinr::read.fasta(seqRNA.loc, as.string = TRUE)
                SS.seq  <- seqinr::read.fasta(seqSS.loc,  as.string = TRUE)
                if(!all(names(RNA.seq) == names(SS.seq))) stop("Cannot match the names of two files.")
                RNA.row <- sapply(RNA.seq, seqinr::getSequence, as.string = TRUE)
                SS.row  <- sapply(SS.seq,  seqinr::getSequence, as.string = TRUE)
                if(withMFE) {
                        seq.ss  <- mapply(function(RNA, SS) {
                                X1 <- RNA
                                X3 <- as.numeric(substr(SS, nchar(RNA) + 3, nchar(SS) -1))
                                X2 <- substr(SS, 1, nchar(RNA))
                                out <- c(X1, X2, X3)
                        }, RNA.row, SS.row)
                        out <- data.frame(seq.ss, stringsAsFactors = FALSE)
                } else {
                        RNA.Seq <- unlist(RNA.row)
                        SS.Seq  <- unlist(SS.row)
                        out <- data.frame(rbind(RNA.Seq, SS.Seq), stringsAsFactors = FALSE)
                }

        } else {
                readFile <- scan(oneFile.loc, what = character(), sep = "\n", quiet = TRUE)

                # if(length(readFile) > 3) message("If separateFile = FALSE, only one sequence will be read.")
                # if(substr(readFile[[1]], 1, 1) == ">") seqName <- substring(readFile[[1]], 2)

                if(withMFE) {
                        X1 <- readFile[[1]]
                        X3 <- as.numeric(substr(readFile[[2]], nchar(readFile[[1]]) + 3, nchar(readFile[[2]]) -1))
                        X2 <- substr(readFile[[2]], 1, nchar(readFile[[1]]))
                        out <- rbind(X1, X2, X3)
                        out <- data.frame(sequence = out, stringsAsFactors = FALSE)
                } else out <- data.frame(sequence = readFile, stringsAsFactors = FALSE)
        }
        out
}

###### make_frequencies ######

#' Make the frequencies file for new classifier construction
#' @description This function is used to calculate the frequencies of lncRNAs, CDs, and
#' secondary structure sequences. The frequencies file can be used to build the classifier
#' using function \code{\link{extract_features}}. Functions \code{make_frequencies} and
#' \code{extract_features} are useful when users are trying
#' to build their own model.
#'
#' NOTE: Function \code{make_frequencies} makes the frequencies file
#' for building the classifiers of LncFinder method. If users need to calculate Logarithm-Distance,
#' Euclidean-Distance, and hexamer score, the frequencies file need to be computed using function
#' \code{\link{make_referFreq}}.
#'
#' @param cds.seq Coding sequences (mRNA without UTRs). Can be a FASTA file loaded
#' by \code{\link[seqinr]{seqinr-package}} or secondary structure
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
#' \code{\link[seqinr]{seqinr-package}} or secondary structure
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
#' @author HAN Siyu
#'
#' @details This function is used to make frequencies file for LncFinder method. This file is needed
#' when users are trying to build their own model.
#'
#' In order to achieve high accuracy, mRNA should not be regarded as CDs and assigned
#' to parameter \code{cds.seq}. However, CDs of some species may be insufficient
#' for calculating frequencies, and mRNAs can be regarded as CDs with parameter
#' \code{check.cds = TRUE}. In this case, hexamer frequencies will be calculated
#' on ORF region.
#'
#' Considering that it is time consuming to obtain secondary structure sequences,
#' users can only provide nucleotide sequences and build a model without secondary
#' structure features (\code{SS.features = } \code{FALSE}). If users want to build a model
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
#' 3 nt each step.
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
#' @importFrom seqinr s2c
#' @importFrom seqinr count
#'
#' @section References:
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
#'
#' @seealso \code{\link{run_RNAfold}}, \code{\link{read_SS}},
#'          \code{\link{build_model}}, \code{\link{extract_features}},
#'          \code{\link{make_referFreq}}.
#' @examples
#' ### Only for examples:
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
#'                               lncRNA.seq = SS.seq, SS.features = TRUE,
#'                               cds.format = "SS", lnc.format = "SS",
#'                               check.cds = TRUE, ignore.illegal = FALSE)
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
                                 step = 3, freq = F, ignore.illegal = ignore.illegal)
        message("  Completed.", "\n")

        message("  Coding sequences")
        cds.DNA.freq <- get.freq(cds.DNA, alphabet = c("a", "c", "g", "t"), wordsize = 6,
                                 step = 3, freq = F, ignore.illegal = ignore.illegal)
        message("  Completed.", "\n")

        message("+ Frequencies calculation completed. ", Sys.time())

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

#' Extract the Features
#' @description This function can construct the dataset. This function is only used
#' to extract the features, please use function \code{\link{build_model}} to build
#' new models.
#'
#' @param Sequences mRNA sequences or long non-coding sequences. Can be a FASTA
#' file loaded by \code{\link[seqinr]{seqinr-package}} or
#' secondary structure sequences (Dot-Bracket Notation) obtained from function
#' \code{\link{run_RNAfold}}. If \code{Sequences} are secondary structure
#' sequences file, parameter \code{format} should be defined as \code{"SS"}.
#'
#' @param label Optional. String. Indicate the label of the sequences such as
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
#' or \code{"wheat"} to use pre-build frequencies files. Or assign a users' own
#' frequencies file (See function \code{\link{make_frequencies}}).
#'
#' @param parallel.cores Integer. The number of cores for parallel computation.
#' By default the number of cores is \code{2}. Users can set as \code{-1} to run
#' this function with all cores.
#'
#' @return Returns a data.frame. 11 features when \code{SS.features} is \code{FALSE},
#' and 19 features when \code{SS.features} is \code{TRUE}.
#'
#' @author HAN Siyu
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
#' other species, \code{SS.features} as \code{TRUE} may lead to low accuracy.
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
#'    Signal at 1/3 position (\code{Signal.Peak});
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
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
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
#'          \code{\link{make_frequencies}}, \code{\link{run_RNAfold}}, \code{\link{read_SS}}.
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' ### Extract features with pre-build frequencies.file:
#' my_features <- extract_features(Seqs, label = "Class.of.the.Sequences",
#'                                 SS.features = FALSE, format = "DNA",
#'                                 frequencies.file = "mouse",
#'                                 parallel.cores = 2)
#'
#' ### Use your own frequencies file by assign frequencies list to parameter
#' ### "frequencies.file".
#' }
#' @export

extract_features <- function(Sequences, label = NULL, SS.features = FALSE, format = "DNA",
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
                } else stop("- Error: Wrong frequencies.file name.")
        } else  Internal.data = frequencies.file

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")

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
        message("- ORF and LogDist.Seq k = 6")
        Seq.ORF  <- parallel::parSapply(cl, DNA.seq, LogDist.DNA, wordsize = 6, step = 3,
                                        hexamer.lnc = Internal.data$DNA.lnc,
                                        hexamer.cds = Internal.data$DNA.cds)

        message("- EIIP", "\n")
        EIIP.DFT <- parallel::parSapply(cl, DNA.seq, EIIP.DFT)
        parallel::stopCluster(cl)
        features <- data.frame(t(Seq.ORF), t(EIIP.DFT))

        message("+ Feature extraction completed. ", Sys.time(), "\n")
        if(SS.features) features <- cbind(features, SS.group)
        if(!is.null(label)) features <- cbind(label, features)
        features
}

###### build_model ######

#' Build Users' Own Model
#' @description This function is used to build new models with users' own data.
#'
#' @param lncRNA.seq Long non-coding sequences. Can be a FASTA file loaded by
#' \code{\link[seqinr]{seqinr-package}} or secondary structure
#' sequences file (Dot-Bracket Notation) obtained from function
#' \code{\link{run_RNAfold}}. If \code{lncRNA.seq} is secondary structure
#' sequences file, parameter \code{lncRNA.format} should be defined as \code{"SS"}.
#'
#' @param mRNA.seq mRNA sequences. FASTA file loaded by \code{\link[seqinr]{read.fasta}} or
#' secondary structure sequences (Dot-Bracket Notation) obtained from function
#' \code{\link{run_RNAfold}}. If \code{mRNA.seq} is secondary structure sequences
#' file, parameter \code{mRNA.format} should be defined as \code{"SS"}.
#'
#' @param frequencies.file String or a list obtained from function
#' \code{\link{make_frequencies}}. Input species name \code{"human"},
#' \code{"mouse"} or \code{"wheat"} to use pre-build frequencies files. Or assign
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
#' both \code{mRNA.format} and \code{lncRNA.format} are set as \code{"SS"}, can
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
#' @author HAN Siyu
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
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
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
#' ### Build the model with pre-build frequencies.file:
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
        lnc.data <- extract_features(lncRNA.seq, label = "NonCoding", SS.features = SS.features,
                                     format = lncRNA.format, frequencies.file = frequencies.file,
                                     parallel.cores = parallel.cores)

        message("+ Extract features of protein-coding sequences.")
        pct.data <- extract_features(mRNA.seq, label = "Coding", SS.features = SS.features,
                                     format = mRNA.format, frequencies.file = frequencies.file,
                                     parallel.cores = parallel.cores)

        dataset <- rbind(lnc.data, pct.data)

        svm.mod <- svm_tune(dataset, gamma.range = gamma.range, cost.range = cost.range,
                            seed = seed, folds.num = folds.num, return.model = TRUE,
                            parallel.cores = parallel.cores)

        # message("+ Build the model with the whole dataset.")
        # svm.mod  <- e1071::svm(label ~ ., data = dataset, scale = TRUE, probability = TRUE, kernel = "radial",
        #                 gamma = best.parameters[[1]][[1]], cost = best.parameters[[1]][[2]])

        message("+ Process completed.")
        svm.mod
}

###### lnc_finder ######

#' Long Non-coding RNA Identification
#' @description This function is used to predict sequences are non-coding transcripts
#' or protein-coding transcripts.
#'
#' @param Sequences Unevaluated sequences. Can be a FASTA file loaded by
#' \code{\link[seqinr]{seqinr-package}} or secondary structure sequences
#' (Dot-Bracket Notation) obtained from function \code{\link{run_RNAfold}}. If
#' \code{Sequences} is secondary structure sequences file, parameter \code{format}
#' should be defined as \code{"SS"}.
#'
#' @param SS.features Logical. If \code{SS.features = TRUE}, secondary structure
#' features will be used.
#'
#' @param format String. Define the format of the \code{Sequences}. Can be
#' \code{"DNA"} or \code{"SS"}. \code{"DNA"} for DNA sequences and \code{"SS"}
#' for secondary structure sequences.
#'
#' @param frequencies.file String or a list obtained from function
#' \code{\link{make_frequencies}}. Input species name \code{"human"}, \code{"mouse"}
#' or \code{"wheat"} to use pre-build frequencies files. Or assign a users' own
#' frequencies file (See function \code{\link{make_frequencies}}).
#'
#' @param svm.model String or a svm model obtained from function \code{\link{build_model}}
#' or \code{\link{svm_tune}}. Input species name \code{"human"}, \code{"mouse"}
#' or \code{"wheat"} to use pre-build models. Or assign a users' own model (See
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
#' @author HAN Siyu
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
#' Siyu Han, Yanchun Liang, Qin Ma, Yangyi Xu, Yu Zhang, Wei Du, Cankun Wang & Ying Li.
#' LncFinder: an integrated platform for long non-coding RNA identification utilizing
#' sequence intrinsic composition, structural information, and physicochemical property.
#' \emph{Briefings in Bioinformatics}, 2018, bby065.
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
#'          \code{\link{extract_features}}, \code{\link{run_RNAfold}}, \code{\link{read_SS}}.
#' @examples
#' \dontrun{
#' data(demo_DNA.seq)
#' Seqs <- demo_DNA.seq
#'
#' ### Input one sequence:
#' OneSeq <- Seqs[1]
#' result_1 <- lnc_finder(OneSeq, SS.features = FALSE, format = "DNA",
#'                        frequencies.file = "human", svm.model = "human",
#'                        parallel.cores = 2)
#'
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

        if(class(svm.model)[1] == "character") {
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

        message("+ Predicting...", "\n")
        pred.res <- stats::predict(svm.mod, seq.data, probability = TRUE)
        results  <- data.frame(Pred = pred.res)
        prob.res <- data.frame(Result = results, Coding.Potential = (attr(pred.res, "probabilities")[,2]))
        output   <- cbind(prob.res, seq.data)
        message("+ Process completed.", "\n")
        output
}

