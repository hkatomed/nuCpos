test_predNuCposActLikePredNuPoP <- function(){
    predNuCposActLikePredNuPoP(system.file("extdata", "test.seq", 
        package="nuCpos"), species="mm", smoothHBA=FALSE, 
        std=FALSE)
    results <- read.table(file = "test.seq_Prediction4.txt", skip = 1)
    expect_equal(results$V1[101], 101)
    expect_equal(results$V2[101], 0.001)
    expect_equal(results$V3[101], 0.742)
    expect_equal(results$V4[101], 0)
    expect_equal(results$V5[101], -2.303)
}
