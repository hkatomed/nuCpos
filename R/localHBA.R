localHBA <- function(inseq, species = "mm", silent = FALSE){

    if(silent == FALSE) message("species: ", species, "\n")
    if(!is(inseq)[1] == "character"){
        if(is(inseq)[1] == "DNAString"){
            if(requireNamespace("Biostrings", quietly = TRUE)){
                inseq <- as.character(inseq)
                if(silent == FALSE){
                message(
                "The class of inseq was changed from DNAString to character")
                }
            }else{
            message("DNAString cannot be changed to a character string")
            out <- NA; names(out) <- "HBA"; return(out)
            }
        }else{
            message("The class of inseq must be DNAString or character")
            out <- NA; names(out) <- "HBA"; return(out)
        }
    }
    if(nchar(inseq) != 147){
        message("Length of inseq: ", nchar(inseq), "bp\n")
        message("The length of inseq must be 147 bp")
        out <- NA; names(out) <- "HBA"; return(out)
    }

    if(species == "sc"){
        freqL1 <- nature11142_s2.147.freqL
        tranL1 <- nature11142_s2.147.tranL
        TtranL2 <- nature11142_s2.147.tranL2
        TtranL3 <- nature11142_s2.147.tranL3
        TtranL4 <- nature11142_s2.147.tranL4
        TfreqN4SA <- nature11142_s2.147.freqN4SA
        TfreqN4SB <- nature11142_s2.147.freqN4SB
        TfreqN4SC <- nature11142_s2.147.freqN4SC
        TfreqN4SD <- nature11142_s2.147.freqN4SD
        TfreqN4SE <- nature11142_s2.147.freqN4SE
        TfreqN4SF <- nature11142_s2.147.freqN4SF
        TfreqN4SG <- nature11142_s2.147.freqN4SG
        TfreqN4SH <- nature11142_s2.147.freqN4SH
        TfreqN4SI <- nature11142_s2.147.freqN4SI
        TfreqN4SJ <- nature11142_s2.147.freqN4SJ
        TfreqN4SK <- nature11142_s2.147.freqN4SK
        TfreqN4SL <- nature11142_s2.147.freqN4SL
        TfreqN4SM <- nature11142_s2.147.freqN4SM
        TtranN4 <- nature11142_s2.147.tranN4_SMA
    }
    if(species == "sp"){
        freqL1 <- sd01.147.freqL
        tranL1 <- sd01.147.tranL
        TtranL2 <- sd01.147.tranL2
        TtranL3 <- sd01.147.tranL3
        TtranL4 <- sd01.147.tranL4
        TfreqN4SA <- sd01.147.freqN4SA
        TfreqN4SB <- sd01.147.freqN4SB
        TfreqN4SC <- sd01.147.freqN4SC
        TfreqN4SD <- sd01.147.freqN4SD
        TfreqN4SE <- sd01.147.freqN4SE
        TfreqN4SF <- sd01.147.freqN4SF
        TfreqN4SG <- sd01.147.freqN4SG
        TfreqN4SH <- sd01.147.freqN4SH
        TfreqN4SI <- sd01.147.freqN4SI
        TfreqN4SJ <- sd01.147.freqN4SJ
        TfreqN4SK <- sd01.147.freqN4SK
        TfreqN4SL <- sd01.147.freqN4SL
        TfreqN4SM <- sd01.147.freqN4SM
        TtranN4 <- sd01.147.tranN4_SMA
    }
    if(species == "mm"){            
        freqL1 <- chem.mm9.freqL
        tranL1 <- chem.mm9.tranL
        TtranL2 <- chem.mm9.tranL2
        TtranL3 <- chem.mm9.tranL3
        TtranL4 <- chem.mm9.tranL4
        TfreqN4SA <- chem.mm9.freqN4SA
        TfreqN4SB <- chem.mm9.freqN4SB
        TfreqN4SC <- chem.mm9.freqN4SC
        TfreqN4SD <- chem.mm9.freqN4SD
        TfreqN4SE <- chem.mm9.freqN4SE
        TfreqN4SF <- chem.mm9.freqN4SF
        TfreqN4SG <- chem.mm9.freqN4SG
        TfreqN4SH <- chem.mm9.freqN4SH
        TfreqN4SI <- chem.mm9.freqN4SI
        TfreqN4SJ <- chem.mm9.freqN4SJ
        TfreqN4SK <- chem.mm9.freqN4SK
        TfreqN4SL <- chem.mm9.freqN4SL
        TfreqN4SM <- chem.mm9.freqN4SM
        TtranN4 <- chem.mm9.tranN4_SMA
    }
    outlist <- .Fortran("localHBA_3", 
            inseq = inseq, logascSA = numeric(length=1), 
            logascSB = numeric(length=1), logascSC = numeric(length=1),
            logascSD = numeric(length=1), logascSE = numeric(length=1),
            logascSF = numeric(length=1), logascSG = numeric(length=1),
            logascSH = numeric(length=1), 
            logascSI = numeric(length=1), logascSJ = numeric(length=1),
            logascSK = numeric(length=1), logascSL = numeric(length=1),
            logascSM = numeric(length=1), freqL1 = freqL1, 
            tranL1 = tranL1, TtranL2 = TtranL2, TtranL3 = TtranL3, 
            TtranL4 = TtranL4, TfreqN4SA = TfreqN4SA, TfreqN4SB = TfreqN4SB, 
            TfreqN4SC = TfreqN4SC, TfreqN4SD = TfreqN4SD, 
            TfreqN4SE = TfreqN4SE, 
            TfreqN4SF = TfreqN4SF, TfreqN4SG = TfreqN4SG, 
            TfreqN4SH = TfreqN4SH, TfreqN4SI = TfreqN4SI, 
            TfreqN4SJ = TfreqN4SJ, TfreqN4SK = TfreqN4SK, 
            TfreqN4SL = TfreqN4SL, 
            TfreqN4SM = TfreqN4SM,
            TtranN4 = TtranN4, PACKAGE = "nuCpos")[2:14]

    # N を含んでいるときはすべてを NA にする。
    if(outlist[[1]] == 0){
        outlist[seq_len(13)] <- as.numeric(NA)
    }
    out <- unlist(outlist)
    out <- out[c(1,8,2,9,3,10,4,11,5,12,6,13,7)]
    names(out) <- c("lHBA_A", "lHBA_B", "lHBA_C", "lHBA_D", 
            "lHBA_E", "lHBA_F", "lHBA_G", 
            "lHBA_H", "lHBA_I", "lHBA_J", "lHBA_K", "lHBA_L", "lHBA_M")
    return(out)
}
