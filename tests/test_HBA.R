test_HBA <- function(){
    load(system.file("extdata","inseq.RData",package="nuCpos"))
    load(system.file("extdata","INSEQ_DNAString.RData",package="nuCpos"))
    inseq.N <- gsub(pattern = "A", replacement = "N", inseq)
    mm.HBA <- -5.108546
    sc.HBA <- -2.460025
    sp.HBA <- -2.627370
    expect_equal(unname(HBA(inseq, species = "mm", silent = TRUE)), 
        mm.HBA, tolerance = 1.0e-6)
    expect_equal(unname(HBA(inseq, species = "sc", silent = TRUE)), 
        sc.HBA, tolerance = 1.0e-6)
    expect_equal(unname(HBA(inseq, species = "sp", silent = TRUE)), 
        sp.HBA, tolerance = 1.0e-6)
    expect_equal(unname(HBA(INSEQ, species = "mm", silent = TRUE)), 
        mm.HBA, tolerance = 1.0e-6)
    expect_equal(unname(HBA(INSEQ, species = "sc", silent = TRUE)), 
        sc.HBA, tolerance = 1.0e-6)
    expect_equal(unname(HBA(INSEQ, species = "sp", silent = TRUE)), 
        sp.HBA, tolerance = 1.0e-6)
    expect_true(is.na(HBA("AAA", species = "mm", silent = TRUE)))
    expect_true(is.na(HBA(123, species = "mm", silent = TRUE)))
    expect_true(is.na(HBA(inseq.N, species = "mm", silent = TRUE)))
}
