predNuCposActLikePredNuPoP <- function (file, species = "mm", smoothHBA = FALSE,
        std = FALSE) {
    file = as.character(file)
    n = nchar(file)
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
        results = .Fortran("nuCpos_2", n, file, freqL, tranL, 
            tranL2, tranL3, tranL4, freqN4, tranN4, maxlen =  as.integer(500), 
            Pd, std = STD, ind = as.integer(0), PACKAGE = "nuCpos")
        ind = results$ind
    }
    if(smoothHBA == TRUE){
        results = .Fortran("nuCpos_1", n, file, freqL, tranL, 
            tranL2, tranL3, tranL4, freqN4, tranN4, maxlen =  as.integer(500), 
            Pd, std = STD, ind = as.integer(0), PACKAGE = "nuCpos")
        ind = results$ind
    }

    if(ind==1){
        warning("The input file specified by the file argument 
        does not exist in the working directory")
    } else if(ind==2){
        warning("The input file must be in single FASTA format.")
        warning("Only the four DNA bases (A/C/G/T) can be accepted.")
        warning("The length of each line must be the same (e.g. 60 bp).")
    } else {
        wd <- getwd()
        filename <- paste(file, "_Prediction4.txt", sep = "")
        message(filename, "was created in the working directory.")
    }
}
