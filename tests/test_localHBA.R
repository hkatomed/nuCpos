test_localHBA <- function(){
    load(system.file("extdata","inseq.RData",package="nuCpos"))
    load(system.file("extdata","INSEQ_DNAString.RData",package="nuCpos"))
    inseq.N <- gsub(pattern = "A", replacement = "N", inseq)
    mm.lHBA <- c(-1.26144039, -1.60878614, 0.04168163, 
        0.67028283, 0.64609413, -2.04965343, -2.87359702, 
        -0.23010702, -0.45807823, -0.46043330, -0.45175477, 
        0.02487367, -0.30991794)
    sc.lHBA <- c(-1.56140949, -1.62502354, 0.48885990, 
        2.37615568, 2.90458625, -1.35195919, -3.13228907, 
        -0.32208031, 0.27650871, 0.01922002, 0.49787625, 
        -0.17151500, -1.27186158)
    sp.lHBA <- c(-1.566757163, -2.249651890, 1.188983606,
        1.808008192, 2.304915648, -0.290338951, -1.741081053, 
        -0.093952092, 0.119058916, -1.335654721, -0.001721381, 
        0.244317796, -0.968842314)
    expect_equal(unname(localHBA(inseq, species = "mm", silent = TRUE)), 
        mm.lHBA, tolerance = 1.0e-6)
    expect_equal(unname(localHBA(inseq, species = "sc", silent = TRUE)), 
        sc.lHBA, tolerance = 1.0e-6)
    expect_equal(unname(localHBA(inseq, species = "sp", silent = TRUE)), 
        sp.lHBA, tolerance = 1.0e-6)
    expect_equal(unname(localHBA(INSEQ, species = "mm", silent = TRUE)), 
        mm.lHBA, tolerance = 1.0e-6)
    expect_equal(unname(localHBA(INSEQ, species = "sc", silent = TRUE)), 
        sc.lHBA, tolerance = 1.0e-6)
    expect_equal(unname(localHBA(INSEQ, species = "sp", silent = TRUE)), 
        sp.lHBA, tolerance = 1.0e-6)
    expect_true(is.na(localHBA("AAA", species = "mm", silent = TRUE)))
    expect_true(is.na(localHBA(123, species = "mm", silent = TRUE)))
    expect_true(is.na(localHBA(inseq.N, species = "mm", silent = TRUE)))
}
