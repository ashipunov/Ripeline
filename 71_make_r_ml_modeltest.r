## this script compares different nucleotide substitution models
## this is very separate script which might not be required to run each time

library(phangorn)

DATE <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp

conc <- as.phyDat(read.dna("40_concatenated/semistrict.fasta", format="fasta"))
LAB <- sub("__.*$", "", labels(conc))

treesp <- read.table("_kubricks_treesp.txt", sep="\t", h=TRUE, as.is=TRUE) # read treesp to understand outlliers
outliers <- treesp$SPECIES.NEW[treesp$TYPE == "outlier" & treesp$USE == 1]
if (length(outliers) > 0) {
 outliers <- gsub(" ", "_", outliers)
 EXC <- labels(conc)[LAB %in% outliers]
 conc <- conc[!labels(conc) %in% EXC, ] # remove outliers, if any
}

models <- modelTest(conc,
 model=c("JC", "HKY", "GTR"), # change the model set if needed
 multicore=TRUE) # multicore=TRUE might not work under Windows, remove if needed

models[order(models$AIC), ] # outputs the table with models arranged by AIC, change this order if needed
