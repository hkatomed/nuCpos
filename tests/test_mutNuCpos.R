test_mutNuCpos <- function(){
    TALS <- paste(scan(file = system.file("extdata", "TALS.fasta", 
        package="nuCpos"), what = character(), skip = 1), sep = "", 
        collapse = "")
    TTAGGGx12 <- paste(scan(file = system.file("extdata", 
        "TTAGGGx12.fasta", package="nuCpos"), what = character(), 
        skip = 1), sep = "", collapse = "")
    results <- mutNuCpos(TALS, site = 1464, ins= TTAGGGx12, species="sc", 
        prob.dyad = TRUE, smoothHBA=TRUE, plot.window = 601, 
        ylim.HBA = c(-11, 0), 
        annotation = data.frame(name = "alpha2", 
        color = "purple", left = 1534, right = 1559), full = TRUE)
    expect_equal(results$results$pstart[101], 0.0005225962, tolerance = 1.0e-8)
    expect_equal(results$results$nucoccup[101], 0.9885672, tolerance = 1.0e-6)
    expect_equal(results$results$viterbi[101], 1)
    expect_equal(results$results$affinity[101], -0.89849, tolerance = 1.0e-5)
    expect_equal(results$results$pos[101], -3673)
}
