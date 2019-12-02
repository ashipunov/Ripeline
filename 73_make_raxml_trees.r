## this script takes the semistrict data and applies RAxML utility to make maximal likelihood tree
## requires RAxML installation an dthe presence of working "raxmlHPC" binary

library(shipunov)
library(ips) # this package contains raxml() function to call RAxML from within R

DATE <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp

treesp <- read.table("_kubricks_treesp.txt", sep="\t", h=TRUE, as.is=TRUE) # is needed to understand outgroups and outliers
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

setwd("70_raxml_working") # go to RAxML working directory
tr <- raxml(conc, m="GTRGAMMAI", # change the model if needed
 f="a", # change RAxML algorithm if needed
 N=100, # change bootstrap options if needed
 p=1234, x=1234, exec="raxmlHPC", # change name of RAxML binary; on Windows, you will _need_ to change it
 outgroup=OUT)
setwd("..") # go out of RAxML working directory
## print "bipartitions" tree to PDF
pdf(paste0("99_trees/", DATE, "_semistrict_raxml_kubricks.pdf"), height=8, width=12) # change PDF size if needed
oldpar <- par(mar=rep(0, 4))
plot(tr$bipartitions)
nodelabels(tr$bipartitions$node.label, frame="none", bg="transparent", adj=-0.1)
mtext("semistrict RAxML", font=2, line=-1)
add.scale.bar()
dev.off()
## also save it as Newick
tr$node.label[tr$node.label == "NA"] <- "" # useful for some Newick reading software
write.tree(tr$bipartitions, file=paste0("99_trees/", DATE, "_semistrict_raxml_kubricks.tre"))
