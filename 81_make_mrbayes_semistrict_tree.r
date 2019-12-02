## this script takes semistric data and applies MrBayes utility to build Bayesian tree
## requires MrBayes installation and working "mb-mpi

library(shipunov)
library(ape)

DATE <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp

treesp <- read.table("_kubricks_treesp.txt", sep="\t", h=TRUE, as.is=TRUE) # to understand outgroups and outliers
outgroups <- treesp$SPECIES.NEW[treesp$TYPE == "outgroup" & treesp$USE == 1]
outgroups <- gsub(" ", "_", outgroups)

conc <- read.dna("40_concatenated/semistrict.fasta", format="fasta")
LAB <- sub("__.*$", "", labels(conc))
OUT <- labels(conc)[LAB %in% outgroups]

outliers <- treesp$SPECIES.NEW[treesp$TYPE == "outlier" & treesp$USE == 1]
if (length(outliers) > 0) {
 outliers <- gsub(" ", "_", outliers)
 EXC <- labels(conc)[LAB %in% outliers]
 conc <- conc[!labels(conc) %in% EXC, ] # remove outliers, if any
}

setwd("80_mrbayes_working") # go to MrBayes working directory
tr <- MrBayes(conc, file="semistrict", exec="mb-mpi", # change MrBayes binary if needed, on Windows, you _need_ to change it
 ngen=1e+04, # change MrBayes options if needed
 run=TRUE) # default is not run, just make a NEXUS file for MrBayes
setwd("..")

tr <- tr[[1]] # we need the first tree
tr <- root(tr, outgroup=OUT, resolve.root=TRUE)
tr$node.label <- suppressWarnings(round(as.numeric(tr$node.label)*100)) # warning is OK so it is suppressed
## plot tree into PDF
pdf(paste0("99_trees/", DATE, "_semistrict_mb_kubricks.pdf"), height=8, width=12) # change PDF size if needed
oldpar <- par(mar=rep(0, 4))
plot(tr)
nodelabels(tr$node.label, frame="none", bg="transparent", adj=-0.1)
mtext("semistrict MB, all compatible to 50% majority rule", font=2, line=-1)
tmp <- legend("bottom", plot=FALSE, legend="") # this is how to get rid of overlapped scale bar
add.scale.bar(x=tmp$text$x, y=tmp$text$y) # it is now centered
dev.off()
## also save it into Newick
tr$node.label[tr$node.label == "NA"] <- "" # useful for some Newick reading software
write.tree(tr, file=paste0("99_trees/", DATE, "_semistrict_mb_kubricks.tre"))
