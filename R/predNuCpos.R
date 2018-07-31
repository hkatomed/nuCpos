predNuCpos <- function(file, inseq, species = "mm", 
            smoothHBA = FALSE, std = FALSE, ActLikePredNuPoP = FALSE){

    if(ActLikePredNuPoP == FALSE){
        results <- predNuCposInternal(inseq, species = species, 
                        smoothHBA = smoothHBA, std = std)
        return(results)
    }
    if(ActLikePredNuPoP == TRUE){
        predNuCposActLikePredNuPoP(file, species = species, 
                        smoothHBA = smoothHBA, std = std)
    }
}
