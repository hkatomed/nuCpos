\name{localHBA}
\alias{localHBA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    R function for calculating the histone binding score 
    of a given 147-bp sequence.
}
\description{
    This function invokes a Fortran subroutine to calculate 
    local histone binding score. Nucleosomal and linker 
    models built upon the chemical maps are used for the 
    calculation.}
\usage{localHBA(inseq, species = "mm", silent = FALSE)}

\arguments{
    \item{inseq}{a character or DNAString object. 
        The length of the character string must be 147 bp.}
    \item{species}{a character = mm, sc or sp; 
        "mm" for mouse, "sc" for \emph{S. cerevisiae} and 
        "sp" for \emph{S. pombe}.}
    \item{silent}{a logical value indicating whether 
        messages are printed in the console.}
}
\value{\code{localHBA} outputs a numeric vector of length 13: 
    local histone binding affinity scores for 
    specific regions in a nucleosome.
}

\examples{
load(system.file("extdata","inseq.RData",package="nuCpos"))
localHBA(inseq, species = "sc")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

