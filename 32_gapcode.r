## this script appires shipunov::Gap.code() to each trimmed alignmenmt

library(shipunov)

## empty list with names-fragments
dna <- structure(vector("list", 2), names=c("abcd", "efgh"))

## go to directory with trimmed alignments
setwd("31_alignments_trimmed")

for (f in names(dna)) {
cat(f, "\n")
dna[[f]] <- Read.fasta(paste0(f, ".fasta"))
gc <- Gap.code(dna[[f]]$sequence) # Gap.code() takes character vector where each element is a sequence string
dna[[f]]$sequence <- apply(cbind(dna[[f]]$sequence, gc), 1, paste, collapse="") # Gap.code() outputs character _matrix_ where each column is a gapcoded posititon
Write.fasta(dna[[f]], file=paste0("../32_alignments_trimmed_gapcoded/", f, ".fasta"))
}
