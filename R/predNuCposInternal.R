predNuCposInternal <- function(inseq, species = "mm", 
            smoothHBA = FALSE, std = FALSE){
    
    message("species: ", species)
    if(!is(inseq)[1] == "character"){
        if(is(inseq)[1] == "DNAString"){
            if(requireNamespace("Biostrings", quietly = TRUE)){
                inseq <- as.character(inseq)
                message(paste("The class of inseq was changed ", 
                    "from DNAString to character\n", sep = ""))
            }else{
                stop("DNAString cannot be changed to a character string")
            }
        }else{
            stop("The class of inseq must be DNAString or character")
        }
    }

    inseq <- toupper(inseq)

    if(nchar(inseq) <= 1000){
        message("Length of inseq: ", nchar(inseq), "bp")
        stop("The length of inseq must not be less than 1,000 bp")
    }
    
    inseqlen <- nchar(inseq)

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

    if(smoothHBA == FALSE){
        results = .Fortran("nuCpos2_2", inseq, inseqlen, freqL, tranL, 
            tranL2, tranL3, tranL4, freqN4, tranN4, maxlen =  as.integer(500), 
            Pd, std = STD, 
            pstart = numeric(length=inseqlen),
            nucoccup = numeric(length=inseqlen),
            viterbi = numeric(length=inseqlen),
            affinity = numeric(length=inseqlen),
            PACKAGE = "nuCpos")[13:16]
    }
    if(smoothHBA == TRUE){
        results = .Fortran("nuCpos2_1", inseq, inseqlen, freqL, tranL, 
            tranL2, tranL3, tranL4, freqN4, tranN4, maxlen =  as.integer(500), 
            Pd, std = STD, 
            pstart = numeric(length=inseqlen),
            nucoccup = numeric(length=inseqlen),
            viterbi = numeric(length=inseqlen),
            affinity = numeric(length=inseqlen),
            PACKAGE = "nuCpos")[13:16]
    }
    results <- data.frame(results)
    results$pos <- seq_len(inseqlen)
    results <- results[,c(5, 1, 2, 3, 4)]
    results$affinity[seq_len(73)] <- as.numeric(NA)
    results$affinity[(inseqlen-72):inseqlen] <- as.numeric(NA)
    return(results)
}
