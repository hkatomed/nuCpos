mutNuCpos <- function(wtseq, site = 1, ins = "", del = 0, species = "mm", 
            smoothHBA = FALSE, std = FALSE, plot.window = 501, 
            prob.dyad = FALSE, show.viterbi = FALSE, occup.window = 200, 
            show.occup.window = FALSE, ymax.prob = 1.1, ymax.occup = 1.1, 
            ylim.HBA = c(-15, 5), annotation = data.frame(name = "", 
            color = "", left = 0, right = 0)[0,], full = FALSE){

    if(plot.window%%2 == 0) stop("\"plot.window\" must be an odd number.")
    if(occup.window%%2 == 1) stop("\"occup.window\" must be an even number.")
    ylim.change <- FALSE
    if(del < 0) stop("\"del\" must be >0.")
    if(nchar(ins) > plot.window){
        stop("Length of \"ins\" must not be longer than \"plot.window\".")
    }
    
    message("species: ", species)
    if(!is(wtseq)[1] == "character"){
        if(is(wtseq)[1] == "DNAString"){
            if(requireNamespace("Biostrings", quietly = TRUE)){
                wtseq <- as.character(wtseq)
                message(paste("The class of wtseq was changed ", 
                    "from DNAString to character\n", sep = ""))
            }else{
                stop("DNAString cannot be changed to a character string")
            }
        }else{
            stop("The class of wtseq must be DNAString or character")
        }
    }
    if(!is(ins)[1] == "character"){
        if(is(ins)[1] == "DNAString"){
            if(requireNamespace("Biostrings", quietly = TRUE)){
                ins <- as.character(ins)
                message(paste("The class of wtseq was changed ", 
                    "from DNAString to character\n", sep = ""))
            }else{
                stop("DNAString cannot be changed to a character string")
            }
        }else{
            stop("The class of ins must be DNAString or character")
        }
    }

    wtseq <- toupper(wtseq)
    ins <- toupper(ins)

    if(nchar(wtseq) <= 1000){
        message("Length of wtseq: ", nchar(wtseq), "bp")
        stop("The length of wtseq must not be less than 1,000 bp")
    }else{
        message("Length of wtseq: ", nchar(wtseq), "bp")
        leftseq <- substring(wtseq, first = 1, last = site - 1)
        rightseq <- substring(wtseq, first = site + del, last = nchar(wtseq))
        mtseq <- paste(leftseq, ins, rightseq, sep = "", collapse = "")
        message("Length of mtseq: ", nchar(mtseq), "bp")
        leftflank <- substring(wtseq, first = site - 1 - 9, last = site - 1)
        rightflank <- substring(wtseq, first = site + del, 
            last = site + del + 9)
        rightflank.wt <- substring(wtseq, first = site, last = site + 9)
        if(ins != "") insseq <- ins
        if(ins == "") insseq <- "*"
        message("Wild type: ", site - 1 - 9, "-", 
            leftflank, rightflank.wt, sep = "")
        message("Mutant:    ", site - 1 - 9, "-", 
            leftflank, "[", insseq, "]", rightflank, sep = "")
    }
    wtseqx5 <- paste(wtseq, wtseq, wtseq, wtseq, wtseq, sep = "", collapse = "")
    mtseqx5 <- paste(mtseq, mtseq, mtseq, mtseq, mtseq, sep = "", collapse = "")

    wtlen <- nchar(wtseqx5)
    mtlen <- nchar(mtseqx5)

    if(std == TRUE) STD <- 1
    if(std == FALSE) STD <- 0

    if(species == "sc"){
        freqL <- nature11142_s2.147.freqL
        tranL <- nature11142_s2.147.tranL
        tranL2 <- nature11142_s2.147.tranL2
        tranL3 <- nature11142_s2.147.tranL3
        tranL4 <- nature11142_s2.147.tranL4
        freqN4 <- nature11142_s2.147.freqN4SA
        tranN4 <- nature11142_s2.147.tranN4_SMA
        Pd <- nature11142_s2.linker.147.prob_SMA
    }
    if(species == "sp"){
        freqL <- sd01.147.freqL
        tranL <- sd01.147.tranL
        tranL2 <- sd01.147.tranL2
        tranL3 <- sd01.147.tranL3
        tranL4 <- sd01.147.tranL4
        freqN4 <- sd01.147.freqN4SA
        tranN4 <- sd01.147.tranN4_SMA
        Pd <- sd01.linker.147.prob_SMA
    }
    if(species == "mm"){            
        freqL <- chem.mm9.freqL
        tranL <- chem.mm9.tranL
        tranL2 <- chem.mm9.tranL2
        tranL3 <- chem.mm9.tranL3
        tranL4 <- chem.mm9.tranL4
        freqN4 <- chem.mm9.freqN4SA
        tranN4 <- chem.mm9.tranN4_SMA
        Pd <- chem.mm9.LinkerDNA.prob_SMA
    }

    ## wtseq
    if(smoothHBA == FALSE){
        results = .Fortran("nuCpos2_2", wtseqx5, wtlen, freqL, tranL, 
            tranL2, tranL3, tranL4, freqN4, tranN4, maxlen =  as.integer(500), 
            Pd, std = STD, 
            pstart = numeric(length=wtlen),
            nucoccup = numeric(length=wtlen),
            viterbi = numeric(length=wtlen),
            affinity = numeric(length=wtlen),
            PACKAGE = "nuCpos")[13:16]
    }
    if(smoothHBA == TRUE){
        results = .Fortran("nuCpos2_1", wtseqx5, wtlen, freqL, tranL, 
            tranL2, tranL3, tranL4, freqN4, tranN4, maxlen =  as.integer(500), 
            Pd, std = STD, 
            pstart = numeric(length=wtlen),
            nucoccup = numeric(length=wtlen),
            viterbi = numeric(length=wtlen),
            affinity = numeric(length=wtlen),
            PACKAGE = "nuCpos")[13:16]
    }
    wt.results <- data.frame(results)
    wt.results$affinity[seq_len(73)] <- as.numeric(NA)
    wt.results$affinity[(wtlen-72):wtlen] <- as.numeric(NA)

    ## mtseq
    if(smoothHBA == FALSE){
        results = .Fortran("nuCpos2_2", mtseqx5, mtlen, freqL, tranL, 
            tranL2, tranL3, tranL4, freqN4, tranN4, maxlen =  as.integer(500), 
            Pd, std = STD, 
            pstart = numeric(length=mtlen),
            nucoccup = numeric(length=mtlen),
            viterbi = numeric(length=mtlen),
            affinity = numeric(length=mtlen),
            PACKAGE = "nuCpos")[13:16]
    }
    if(smoothHBA == TRUE){
        results = .Fortran("nuCpos2_1", mtseqx5, mtlen, freqL, tranL, 
            tranL2, tranL3, tranL4, freqN4, tranN4, maxlen =  as.integer(500), 
            Pd, std = STD, 
            pstart = numeric(length=mtlen),
            nucoccup = numeric(length=mtlen),
            viterbi = numeric(length=mtlen),
            affinity = numeric(length=mtlen),
            PACKAGE = "nuCpos")[13:16]
    }
    mt.results <- data.frame(results)
    mt.results$affinity[seq_len(73)] <- as.numeric(NA)
    mt.results$affinity[(mtlen-72):mtlen] <- as.numeric(NA)


    old.par <- par(no.readonly=TRUE)
    if(show.viterbi == TRUE)    par(mfrow = c(4,2))
    if(show.viterbi == FALSE)    par(mfrow = c(3,2))
    
    wt.site <- nchar(wtseq) * 2 + site
    wt.from <- wt.site - round(plot.window/2)
    wt.to <- wt.site + round(plot.window/2)
    wt.index <- -round(plot.window/2):round(plot.window/2)
    mt.site <- nchar(mtseq) * 2 + site
    mt.from <- mt.site - round(plot.window/2)
    mt.to <- mt.site + round(plot.window/2) + nchar(ins) - del
    mt.index <- -round(plot.window/2):(round(plot.window/2) + nchar(ins) - del)

    # Difference in occupancy
    wt.occu.left.from <- wt.site - round(occup.window/2)
    wt.occu.left.to <- wt.site - 1
    wt.occu.right.from <- wt.site + 1 + del
    wt.occu.right.to <- wt.site + 1 + del + round(occup.window/2)

    mt.occu.left.from <- mt.site - round(occup.window/2)
    mt.occu.left.to <- mt.site - 1
    mt.occu.right.from <- mt.site + 1 + nchar(ins)
    mt.occu.right.to <- mt.site + 1 + nchar(ins) + round(occup.window/2)

    wt.occu.left <- wt.results$nucoccup[wt.occu.left.from:wt.occu.left.to]
    wt.occu.right <- wt.results$nucoccup[wt.occu.right.from:wt.occu.right.to]
    mt.occu.left <- mt.results$nucoccup[mt.occu.left.from:mt.occu.left.to]
    mt.occu.right <- mt.results$nucoccup[mt.occu.right.from:mt.occu.right.to]
    wt.occu <- c(wt.occu.left, wt.occu.right)
    mt.occu <- c(mt.occu.left, mt.occu.right)
    diff.length <- round(occup.window/2) + round(occup.window/2)
    message("occup.window: ", diff.length)
    diff.occu <- sum(abs(wt.occu - mt.occu))
    names(diff.occu) <- paste("Difference in occupancy (", occup.window, 
                " bp)", sep = "")

    xmin <- -round(plot.window/2)
    xmax <- max(round(plot.window/2), round(plot.window/2) + nchar(ins) - del)
    xlim <- c(xmin, xmax)

    if(annotation$left[1] == 0 && annotation$right[1] == 0){
        draw.annotation <- FALSE
    }else{
        draw.annotation <- TRUE
        annotation$name <- as.character(annotation$name)
        annotation$color <- as.character(annotation$color)
        message("----------")
        message("annotation\n", "name color left right")
        for(i in seq_len(nrow(annotation))){
            message(paste("   ", annotation[i,], sep = "", collapse = ""))
        }
        message("----------")
    }
    
    xlab <- "Distance from the target site (bp)"

    # Probability
    ylim <- c(-0.1, ymax.prob)
    if(prob.dyad == FALSE){
        plot(x = wt.index, y = wt.results$pstart[wt.from:wt.to], 
            type = "h", col = "blue", lwd = 1, xlab = xlab, 
            ylab = "Probability", xlim = xlim, ylim = ylim, yaxt = "n")
    }
    if(prob.dyad == TRUE){
        plot(x = wt.index, y = wt.results$pstart[(wt.from - 73):(wt.to - 73)], 
            type = "h", col = "blue", lwd = 1, xlab = xlab, 
            ylab = "Probability", xlim = xlim, ylim = ylim, yaxt = "n")
    }
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), las = 2, 
        labels = c("0", "", "0.5", "", "1"))
    if(del > 0){
        rect(xleft = 0, ybottom = -0.08, xright = del, ytop = -0.05, 
            col = "black", border = "black")
    }
    if(prob.dyad == FALSE)    title(main = "5'-end probability: Wild type")
    if(prob.dyad == TRUE)    title(main = "Dyad probability: Wild type")
    
    if(draw.annotation == TRUE){
        for(i in seq_len(nrow(annotation))){
            xleft <- annotation$left[i] - site
            xright <- annotation$right[i] - site
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = -0.08, xright = xright, ytop = -0.05, 
                col = color, border = color)
        }
    }

    if(prob.dyad == FALSE){
        plot(x = mt.index, y = mt.results$pstart[mt.from:mt.to], 
            type = "h", col = "blue", lwd = 1, xlab = xlab, 
            ylab = "Probability", xlim = xlim, ylim = ylim, yaxt = "n")
    }
    if(prob.dyad == TRUE){
        plot(x = mt.index, y = mt.results$pstart[(mt.from - 73):(mt.to - 73)], 
            type = "h", col = "blue", lwd = 1, xlab = xlab, 
            ylab = "Probability", xlim = xlim, ylim = ylim, yaxt = "n")
    }
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), las = 2, 
        labels = c("0", "", "0.5", "", "1"))
    if(nchar(ins) > 0){
        rect(xleft = 0, ybottom = -0.08, xright = nchar(ins), ytop = -0.05, 
            col = "red", border = "red")
    }
    if(prob.dyad == FALSE)    title(main = "5'-end probability: Mutant")
    if(prob.dyad == TRUE)    title(main = "Dyad probability: Mutant")
    
    if(draw.annotation == TRUE){
        for(i in seq_len(nrow(annotation))){
            if(annotation$left[i] > site){
                xleft <- annotation$left[i] - site + nchar(ins) - del
                xright <- annotation$right[i] - site + nchar(ins) - del
            }else{
                xleft <- annotation$left[i] - site
                xright <- annotation$right[i] - site
            }
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = -0.08, xright = xright, ytop = -0.05, 
                col = color, border = color)
        }
    }
    
    
    # Occupancy
    ylim <- c(-0.1, ymax.occup)
    plot(x = wt.index, y = wt.results$nucoccup[wt.from:wt.to], type = "n", 
        yaxt = "n", xlab = xlab, ylab = "Occupancy", xlim = xlim, ylim = ylim)
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), las = 2, 
        labels = c("0", "", "0.5", "", "1"))
    polygon(x = c(wt.index[1], wt.index, wt.index[length(wt.index)]), 
        y = c(0, wt.results$nucoccup[wt.from:wt.to], 0), col = "gray")
    if(del > 0){
        rect(xleft = 0, ybottom = -0.08, xright = del, ytop = -0.05, 
            col = "black", border = "black")
    }
    if(show.occup.window == TRUE){
        polygon(x = c(-round(occup.window/2), -round(occup.window/2):-1, -1), 
            y = c(0, wt.occu.left, 0), density = 20, col = "black")
        polygon(x = c(1 + del, 
                (1 + del):(1 + del + round(occup.window/2)), 
                1 + del + round(occup.window/2)), 
            y = c(0, wt.occu.right, 0), density = 20, col = "black")
    }
    title(main = "Nucleosome occupancy: Wild type")

    if(draw.annotation == TRUE){
        for(i in seq_len(nrow(annotation))){
            xleft <- annotation$left[i] - site
            xright <- annotation$right[i] - site
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = -0.08, xright = xright, ytop = -0.05, 
                col = color, border = color)
        }
    }

    plot(x = mt.index, y = mt.results$nucoccup[mt.from:mt.to], type = "n", 
        yaxt = "n", xlab = xlab, ylab = "Occupancy", xlim = xlim, ylim = ylim)
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), las = 2, 
        labels = c("0", "", "0.5", "", "1"))
    polygon(x = c(mt.index[1], mt.index, mt.index[length(mt.index)]), 
        y = c(0, mt.results$nucoccup[mt.from:mt.to], 0), col = "gray")
    if(nchar(ins) > 0){
        rect(xleft = 0, ybottom = -0.08, xright = nchar(ins), ytop = -0.05, 
            col = "red", border = "red")
    }
    if(show.occup.window == TRUE){
        polygon(x = c(-round(occup.window/2), -round(occup.window/2):-1, -1), 
            y = c(0, mt.occu.left, 0), density = 20, col = "black")
        polygon(x = c(1 + nchar(ins), 
                (1 + nchar(ins)):(1 + nchar(ins) + round(occup.window/2)), 
                1 + nchar(ins) + round(occup.window/2)), 
            y = c(0, mt.occu.right, 0), density = 20, col = "black")
    }
    title(main = "Nucleosome occupancy: Mutant")

    if(draw.annotation == TRUE){
        for(i in seq_len(nrow(annotation))){
            if(annotation$left[i] > site){
                xleft <- annotation$left[i] - site + nchar(ins) - del
                xright <- annotation$right[i] - site + nchar(ins) - del
            }else{
                xleft <- annotation$left[i] - site
                xright <- annotation$right[i] - site
            }
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = -0.08, xright = xright, ytop = -0.05, 
                col = color, border = color)
        }
    }


    # Histone Binding Affinity
    ymin <- min(wt.results$affinity[wt.from:wt.to], 
        mt.results$affinity[mt.from:mt.to])
    ymax <- max(wt.results$affinity[wt.from:wt.to], 
        mt.results$affinity[mt.from:mt.to])

    if(ylim.HBA[1] != 0 || ylim.HBA[2] != 0){
        ymin <- ylim.HBA[1]
        ymax <- ylim.HBA[2]
    }
    ylim <- c(ymin - 2, ymax + 2)

    plot(x = wt.index, y = wt.results$affinity[wt.from:wt.to], type = "l", 
        col = "blue", lwd = 1.5, yaxt = "n", 
        xlab = xlab, ylab = "HBA", 
        xlim = xlim, ylim = ylim)
    axis(side = 2, 
        at = c(-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30), 
        las = 2, labels = c("-30", "-25", "-20", "-15", "-10", "-5", "0", 
                "5", "10", "15", "20", "25", "30"))
    if(del > 0){
        rect(xleft = 0, ybottom = ymin - 1.6, xright = del, ytop = ymin - 1.2, 
            col = "black", border = "black")
    }
    title(main = "Histone Binding Affinity: Wild type")

    if(draw.annotation == TRUE){
        for(i in seq_len(nrow(annotation))){
            xleft <- annotation$left[i] - site
            xright <- annotation$right[i] - site
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = ymin - 1.6, xright = xright, 
                ytop = ymin - 1.2, col = color, border = color)
        }
    }

    plot(x = mt.index, y = mt.results$affinity[mt.from:mt.to], type = "l", 
        col = "blue", lwd = 1.5, yaxt = "n", 
        xlab = xlab, ylab = "HBA", 
        xlim = xlim, ylim = ylim)
    axis(side = 2, 
        at = c(-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30), 
        las = 2, labels = c("-30", "-25", "-20", "-15", "-10", "-5", "0", 
                "5", "10", "15", "20", "25", "30"))
    if(nchar(ins)){
        rect(xleft = 0, ybottom = ymin - 1.6, xright = nchar(ins), 
            ytop = ymin - 1.2, col = "red", border = "red")
    }
    title(main = "Histone Binding Affinity: Mutant")
    
    if(draw.annotation == TRUE){
        for(i in seq_len(nrow(annotation))){
            if(annotation$left[i] > site){
                xleft <- annotation$left[i] - site + nchar(ins) - del
                xright <- annotation$right[i] - site + nchar(ins) - del
            }else{
                xleft <- annotation$left[i] - site
                xright <- annotation$right[i] - site
            }
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = ymin - 1.6, xright = xright, 
                ytop = ymin - 1.2, col = color, border = color)
        }
    }


    # Viterbi path
    if(show.viterbi == TRUE){
    ylim <- c(-0.3, 1.3)
    plot(x = wt.index, y = wt.results$viterbi[wt.from:wt.to], type = "l", 
        yaxt = "n", xlab = xlab, ylab = "State", col = "blue", lwd = 1.5, 
        xlim = xlim, ylim = ylim)
    axis(side = 2, at = c(0, 1), las = 2, 
        labels = c("L", "N"))
    title(main = "Viterbi path: Wild type")
    if(del > 0){
        rect(xleft = 0, ybottom = -0.2, xright = del, ytop = -0.25, 
            col = "black", border = "black")
    }

    if(draw.annotation == TRUE){
        for(i in seq_len(nrow(annotation))){
            xleft <- annotation$left[i] - site
            xright <- annotation$right[i] - site
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = -0.2, xright = xright, 
                ytop = -0.25, col = color, border = color)
        }
    }

    plot(x = mt.index, y = mt.results$viterbi[mt.from:mt.to], type = "l", 
        yaxt = "n", xlab = xlab, ylab = "State", col = "blue", lwd = 1.5, 
        xlim = xlim, ylim = ylim)
    axis(side = 2, at = c(0, 1), las = 2, 
        labels = c("L", "N"))
    title(main = "Viterbi path: Mutant")
    if(nchar(ins) > 0){
        rect(xleft = 0, ybottom = -0.2, xright = nchar(ins), ytop = -0.25, 
            col = "red", border = "red")
    }

    if(draw.annotation == TRUE){
        for(i in seq_len(nrow(annotation))){
            if(annotation$left[i] > site){
                xleft <- annotation$left[i] - site + nchar(ins) - del
                xright <- annotation$right[i] - site + nchar(ins) - del
            }else{
                xleft <- annotation$left[i] - site
                xright <- annotation$right[i] - site
            }
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = -0.2, xright = xright, 
                ytop = -0.25, col = color, border = color)
        }
    }

    } # if(show.viterbi == TRUE)

    par(old.par)


    # output
    wt.HBA <- wt.results$affinity[wt.site]
    mt.HBA.left <- mt.results$affinity[mt.site]
    mt.HBA.center <- mt.results$affinity[mt.site + round(nchar(ins)/2)]
    mt.HBA.right <- mt.results$affinity[mt.site + nchar(ins)]
    names(wt.HBA) <- "wild type HBA at position 0"
    names(mt.HBA.left) <- "mutant HBA at position 0"
    names(mt.HBA.center) <- paste("mutant HBA at position ", 
        round(nchar(ins)/2), sep = "")
    names(mt.HBA.right) <- paste("mutant HBA at position ", 
        nchar(ins), sep = "")
    out1 <- c(diff.occu, wt.HBA, mt.HBA.left, mt.HBA.center, mt.HBA.right)

    if(full == FALSE){
        return(out1)
    }
    if(full == TRUE){
        for(i in seq_len(length(out1))){
            message(names(out1)[i], ": ", round(out1[i], digits = 4), 
                sep = "")
        }
        # cat(names(out1))
        # cat(out1)
        mt.results$pos <- (-(nchar(mtseq)*2-1)):(nchar(mtseq)*3)
        mt.results <- mt.results[,c(5, 1, 2, 3, 4)]
        return(mt.results)
        # arguments <- list(seq = mtseq, site = site, ins = ins, del = del, 
        #     species = species, smoothHBA = smoothHBA, std = std, 
        #     prob.dyad = prob.dyad, annotation = annotation)
        # arguments$annotation$left <- 
        #     arguments$annotation$left + nchar(ins) - del
        # arguments$annotation$right <- 
        #     arguments$annotation$right + nchar(ins) - del
        # out2 <- list(results = mt.results, arguments = arguments, 
        #     out1 = out1)
        # return(out2)
    }
    
    ## merge
    # result <- list(wt.results = wt.results, mt.results = mt.results)
    # return(result)
}
