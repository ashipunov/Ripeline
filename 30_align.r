## this script simply runs the external MUSCLE tool and redirect its output into next directory

## go into directory which contains sets
setwd("20_sets")

## empty list with components names as fragments
dna <- structure(vector("list", 2), names=c("abcd", "efgh"))

## run external tool
for (f in names(dna)) system(paste0("muscle -in ", f, ".fasta", " -out ", "../30_alignments/", f, ".fasta"))

## alternatives:
## MAFFT
## for (f in names(dna)) system(paste0("mafft ", f, ".fasta", " > ", "../30_alignments/", f, "_mafft.fasta"))
## ClustalO
## for (f in names(dna)) system(paste0("clustalo -i ", f, ".fasta", " --force -o ", "../30_alignments/", "_clustalo.fasta", " -v"))
