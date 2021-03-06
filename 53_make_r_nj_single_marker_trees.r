## this script takes trimmed and gaproded _alignments_ and buids NJ tree per each marker

library(ape)
library(shipunov)

DATE <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp

fragments <- c("abcd", "efgh") # fragile, change it if needed

treesp <- read.table("_kubricks_treesp.txt", sep="\t", h=TRUE, as.is=TRUE) # read treesp to determine outgroups and outliers
outgroups <- treesp$SPECIES.NEW[treesp$TYPE == "outgroup" & treesp$USE == 1]
outliers <- treesp$SPECIES.NEW[treesp$TYPE == "outlier" & treesp$USE == 1]

for (f in fragments) {
 cat(f, "\n")
 assign(f, read.dna(paste0("32_alignments_trimmed_gapcoded/", f, ".fasta"), format="fasta"))
 data <- get(f)
 LAB <- sub("__.*$", "", labels(data)) # edit labels to match with treesp
 OUT <- labels(data)[LAB %in% outgroups]
 if (length(outliers) > 0) data <- data[!LAB %in% outliers, ]
 tree <- Root1(njs(dist.dna(data, pairwise.deletion=TRUE)), OUT, resolve.root=TRUE) # Root1() does not fail if unable to root on > 1 outgroup
 pdf(paste0("50_technical_trees/", DATE, "_", f, "_nj_kubricks.pdf"), height=8, width=12) # change PDF size if needed
  par(mar=rep(0, 4))
  plot(tree)
  mtext(f, line=-1, font=2)
 dev.off()
}
