\name{mutNuCpos}
\alias{mutNuCpos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{R function for prediction of nucleosome positioning
on a mutant sequence}
\description{This function plots the results of nucleosome 
    positioning prediction for wild type and mutant sequences
    in a specified window. Nucleosomal and linker models 
    built upon the chemical maps are used for the calculation.
    No file is generated in the current directry. }

\usage{mutNuCpos(wtseq, site = 1, ins = "", del = 0, 
    species = "mm", smoothHBA = FALSE, std = FALSE, 
    plot.window = 501, prob.dyad = FALSE, 
    show.viterbi = FALSE, occup.window = 200, 
    show.occup.window = FALSE, ymax.prob = 1.1, 
    ymax.occup = 1.1, ylim.HBA = c(-15, 5), 
    annotation = data.frame(name = "", color = "", left = 0, 
    right = 0)[0, ], full = FALSE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{wtseq}{a character or DNAString object. 
        The wild-type sequence to be mutated. 
        The string must not contain letters other than 
        "A", "C", "G" or "T."}
    \item{site}{an integer. The site of mutagenesis.}
    \item{ins}{a character or DNAString object. 
        The sequence to be inserted at the "site." 
        The string must not contain letters other than 
        "A", "C", "G" or "T." \code{ins=""}  
        indicates no sequence will be inserted.}
    \item{del}{an integer. 
        The length of the deleted region that starts at 
        the "site." \code{del=0} 
        indicates no sequence will be deleted.}
    \item{species}{a character = mm, sc or sp; 
        "mm" for mouse, "sc" for \emph{S. cerevisiae} and 
        "sp" for \emph{S. pombe}.}
    \item{smoothHBA}{a logical value indicating whether 
        smoothing of histone binding affinity should be 
        applied as in the \code{predNuPoP} function 
        of the parental package \pkg{NuPoP}.}
    \item{std}{a logical value indicating whether 
        standardization should be applied to the 
        histone binding affinity score.}
    \item{plot.window}{an integer. The window to be plotted. 
        This must be an odd number.}
    \item{prob.dyad}{a logical value indicating whether 
        the probability for the predicted dyads is plotted.}
    \item{show.viterbi}{a logical value indicating whether 
        the viterbi path is plotted.}
    \item{occup.window}{an integer. The size of the window 
        for the calculation of occupancy difference. 
        occup.window=200 means that the sum of the absolute 
        occupancy difference for the left-side and 
        right-side 100-bp regions flanking the "site" 
        is caluculated.}
    \item{show.occup.window}{a logical value indicating 
        whether the window for the occupancy difference 
        calculation is shown in the occupancy plots.}
    \item{ymax.prob}{an integer. Specify the upper limit of 
        the y axis of the probability plots.}
    \item{ymax.occup}{an integer. Specify the upper limit of 
        the y axis of the occupancy plots.}
    \item{ylim.HBA}{an integer vector of two values. 
        Specify the lower and upper limits of the y axis of 
        the histone binding affinity plots.}
    \item{annotation}{a data frame. 
        Colored bars can be put under the plots.}
    \item{full}{a logical value indicating whether 
        the calculation results will be returned
        as a data frame object.}
}
\value{
When the \code{full} argument is set as \code{TRUE}, 
the prediction results for the 
mutant sequence will be returned as a data frame object.  
The data frame has five columns as that produced by 
\code{predNuCpos} when its argument \code{ActLikePredNuPoP} 
was set as \code{FALSE}:
    \item{pos}{position in the input DNA sequence}
    \item{pstart}{probability that a nucleosome starts at}
    \item{nucoccup}{nucleosome occupancy score}
    \item{viterbi}{Viterbi path 
        (1 and 0 for the nucleosome and linker states, 
        repsectively)}
    \item{affinity}{histone binding affinity score}

When the \code{full} argument was set as \code{FALSE}, 
this function returns a named numeric vector, in which
the occupancy difference and HBA scores
around the target site are stored. 

When \code{ins=""} and \code{del=0} are applied, 
two wild-type 
sequences are used for the calculation and plotting; 
this yields no difference in the occupancy or HBA.
}

\examples{
# Loading the sequence of TALS, a budding yeast 
# minichromosome.
TALS <- paste(scan(file = system.file("extdata", "TALS.fasta", 
    package="nuCpos"), what = character(), skip = 1), sep = "", 
    collapse = "")
        
# Loading the telomere repeat sequence (hTELx12) 
TTAGGGx12 <- paste(scan(file = system.file("extdata", 
    "TTAGGGx12.fasta", package="nuCpos"), what = character(), 
    skip = 1), sep = "", collapse = "")
mutNuCpos(TALS, site = 1464, ins= TTAGGGx12, species="sc", 
    prob.dyad = TRUE, smoothHBA=TRUE, plot.window = 601, 
    ylim.HBA = c(-11, 0), 
    annotation = data.frame(name = "alpha2", 
        color = "purple", left = 1534, right = 1559))
        
# Loading the telomere repeat isomeric sequence (SI-Ax12) 
TGTAGGx12 <- paste(scan(file = system.file("extdata", 
    "TGTAGGx12.fasta", package="nuCpos"), what = character(), 
    skip = 1), sep = "", collapse = "")
mutNuCpos(TALS, site = 1464, ins= TGTAGGx12, species="sc", 
    prob.dyad = TRUE, smoothHBA=TRUE, plot.window = 601, 
    ylim.HBA = c(-11, 0), 
    annotation = data.frame(name = "alpha2", 
    color = "purple", left = 1534, right = 1559))

# DNA sequences used here are from Ichikawa et al. (2014).
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

