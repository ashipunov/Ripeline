## this script makes k-mer trees directly from DNA database

library(ape)
library(kmer)
library(shipunov)

DATE <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp
dna <- read.table("_kubricks_dna.txt", h=TRUE, sep="\t", as.is=TRUE) # read database

## convert database into FASTA, write it on disc and read again as "DNAbin" object
raw <- dna
raw$ids.tmp <- apply(raw[, c("FRAGMENT", "SPECIES.NEW", "SEQUENCE.ID", "SELECT")], 1, function(.x) (paste(.x, collapse="_")))
Write.fasta(raw[, c("ids.tmp", "SEQUENCE")], file="40_concatenated/raw.fasta")
raw <- read.dna("40_concatenated/raw.fasta", format="fasta")
## now build the k-mer tree from whole data inculding sequences which were de-selected
cat("raw k-mer...\n")
raw.tree <- nj(kdistance(raw, k=7))
raw.tree$edge.length[raw.tree$edge.length < 0] <- 0.00001 # this is to increase visibility, found experimentally
raw.tree$edge.length <- (raw.tree$edge.length + 0.005) * 2 # this is to increase visibility, found experimentally
## plot tree into PDF
pdf(paste0("50_technical_trees/", DATE, "_raw_kmer_nj_kubricks.pdf"), width=12, height=12) # change PDF size if needed
oldpar <- par(mar=rep(0, 4))
plot(raw.tree, no.margin=TRUE)
mtext("k-mer: raw, all", font=2, line=-3)
par(oldpar)
dev.off()

## this is raw data minus de-selected sequences
## make FASTA
rawd <- dna[dna$SELECT == 1, ] # remove deselected
rawd$ids.tmp <- apply(rawd[, c("FRAGMENT", "SPECIES.NEW", "SEQUENCE.ID", "SELECT")], 1, function(.x) (paste(.x, collapse="_")))
Write.fasta(rawd[, c("ids.tmp", "SEQUENCE")], file="40_concatenated/rawd.fasta")
rawd <- read.dna("40_concatenated/rawd.fasta", format="fasta")
## build the tree
cat("rawd k-mer...\n")
rawd.tree <- nj(kdistance(rawd, k=7))
rawd.tree$edge.length[rawd.tree$edge.length < 0] <- 0.00001 # this is to increase visibility, found experimentally
rawd.tree$edge.length <- (rawd.tree$edge.length + 0.005) * 2 # this is to increase visibility, found experimentally
## save tree as PDF
pdf(paste0("50_technical_trees/", DATE, "_rawd_kmer_nj_kubricks.pdf"), width=12, height=12) # change PDF size if needed
oldpar <- par(mar=rep(0, 4))
plot(rawd.tree, no.margin=TRUE)
mtext("k-mer: all minus deselected", font=2, line=-3)
par(oldpar)
dev.off()
