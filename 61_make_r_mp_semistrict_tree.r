## this script takes semistrict data and outputs the maximal parsimony (MP) tree
## it also uses bootstrapped trees to produce majority rule consensus

library(phangorn) # will be used to estimate MP tree, no external tool needed
library(shipunov)

DATE <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp

treesp <- read.table("_kubricks_treesp.txt", sep="\t", h=TRUE, as.is=TRUE) # read treesp to understand outgroups and outlliers
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

conc.phydat <- phyDat(conc) # "phangorn" requires everything in its own format

tr1 <- NJ(dist.hamming(conc.phydat)) # make starting tree
conc.pr <- pratchet(conc.phydat, maxit=100, k=20, start=tr1) # change these parsimony ratchet parameters if needed
if (!is.binary(conc.pr)) conc.pr <- multi2di(conc.pr) # sometimes, trees are not binary
conc.pr <- acctran(conc.pr, conc.phydat)
conc.ptrees <- bootstrap.phyDat(conc.phydat, pratchet, bs=100, multicore=TRUE) # change bootstrap replications if needed
conc.prr <- Root1(conc.pr, OUT, resolve.root=TRUE) # Root1() does not fail if unable to root on > 1 outgroup
## write tree to PDF
pdf(paste0("99_trees/", DATE, "_semistrict_mp_kubricks.pdf"), height=8, width=12)
oldpar <- par(mar=rep(0, 4))
conc.prr.b <- plotBS(conc.prr, conc.ptrees, p=50, type="phylogram")
add.scale.bar()
mtext("semistrict MP", font=2, line=-1)
par(oldpar)
dev.off()
## keep tree as Newick file
conc.prr.b$node.label[conc.prr.b$node.label == "NA"] <- "" # useful for some Newick reading software
write.tree(conc.prr.b, file=paste0("99_trees/", DATE, "_semistrict_mp_kubricks.tre"))

## also save all bootstrapped treea as R object for later use
dput(conc.ptrees, file=paste0("99_trees/", DATE, "_strict_mp_ptrees_kubricks.rd"))
## calculate majority rule consensus
cons.50 <- consensus(conc.ptrees, p=0.5) # change if you need other consensus criteria
cons.50r <- Root1(cons.50, OUT)
## write consensus tree to PDF
pdf(paste0("99_trees/", DATE, "_semistrict_mp_consensus_kubricks.pdf"), height=8, width=12) # change PDF size if needed
oldpar <- par(mar=rep(0, 4))
conc.prr.b <- plot(cons.50r)
add.scale.bar()
mtext("semistrict MP, majority rule consensus", font=2, line=-1)
par(oldpar)
dev.off()
