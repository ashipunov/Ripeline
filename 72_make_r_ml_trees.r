## this script takes semistrict data and estimates maximal likelihood (ML) tree

library(shipunov)
library(phangorn) # no external tool needed

DATE <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp

treesp <- read.table("_kubricks_treesp.txt", sep="\t", h=TRUE, as.is=TRUE) # to find out ougroups and outliers
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

conc <- as.phyDat(conc) # "phangorn" requires its own format
conc.nj  <- NJ(dist.ml(conc)) # starter tree
conc.fit <- pml(conc.nj, data=conc)
conc.fit.opt <- optim.pml(conc.fit,
 model="GTR", optGamma=TRUE, optInv=TRUE, rearrangement="stochastic") # change the model if needed
conc.fit.bs <- bootstrap.pml(conc.fit.opt, bs=100, # change number of boostrap replicates when needed
 optNni=TRUE, multicore=TRUE) # change "multicore" option if it does not work
conc.fit.optr <- root(conc.fit.opt$tree, OUT, resolve.root=TRUE)
## plot tree into PDF
pdf(paste0("99_trees/", DATE, "_semistrict_ml_kubricks.pdf"), height=8, width=12) # change PDF size if needed
oldpar <- par(mar=rep(0, 4))
tr <- plotBS(conc.fit.optr, conc.fit.bs, p=50, type="p")
mtext("semistrict ML", font=2, line=-1)
add.scale.bar()
dev.off()
## and also save it as Newick
tr$node.label[tr$node.label == "NA"] <- "" # useful for some Newick reading software
write.tree(tr, file=paste0("99_trees/", DATE, "_semistrict_ml_kubricks.tre"))
