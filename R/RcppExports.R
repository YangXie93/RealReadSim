# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#'[[Rcpp::export]]
NULL

#'[[Rcpp::export]]
NULL

translateOverlap <- function(c1s, c1e, c2s, c2e, as1, ae1, as2, ae2) {
    .Call(`_RealReadSim_translateOverlap`, c1s, c1e, c2s, c2e, as1, ae1, as2, ae2)
}

evalCoverage <- function(pos, width, sampleID, length, minOverlap, minContigLength, nrOfSamples) {
    .Call(`_RealReadSim_evalCoverage`, pos, width, sampleID, length, minOverlap, minContigLength, nrOfSamples)
}

getIdenticalSeqs <- function(starts1, ends1, starts2, ends2, nm1, nm2, minL = 0L) {
    .Call(`_RealReadSim_getIdenticalSeqs`, starts1, ends1, starts2, ends2, nm1, nm2, minL)
}

getIdenticalSeqsList <- function(names1, starts1, ends1, names2, starts2, ends2, minL = 0L) {
    .Call(`_RealReadSim_getIdenticalSeqsList`, names1, starts1, ends1, names2, starts2, ends2, minL)
}

makeFastaOutput <- function(names, seqs, outFile) {
    .Call(`_RealReadSim_makeFastaOutput`, names, seqs, outFile)
}

sequenceToFastaReads <- function(starts, sequence, meanWidth, newFasta, nameTag) {
    .Call(`_RealReadSim_sequenceToFastaReads`, starts, sequence, meanWidth, newFasta, nameTag)
}

calcMinOverlap <- function(seq, meanWidth) {
    .Call(`_RealReadSim_calcMinOverlap`, seq, meanWidth)
}

subSeqs <- function(seq, starts, ends) {
    .Call(`_RealReadSim_subSeqs`, seq, starts, ends)
}

calcCovVec <- function(readsPerSample, lengths) {
    .Call(`_RealReadSim_calcCovVec`, readsPerSample, lengths)
}

