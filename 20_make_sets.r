## this sctipt takes DNA database and for each fragment, make separate FASTA file
## all "our" sequences used
## only those external source sequences used which (a) belong to species / fragments absent in "our" and (b) longest

library(shipunov)

set <- read.table("_kubricks_dna.txt", sep="\t", h=TRUE, as.is=TRUE)
set <- set[order(set$FRAGMENT, set$SPECIES.NEW), ]
## do not use deselected sequences
set <- set[set$SELECT == 1 | is.na(set$SELECT), ]

## split data frame into list
sets <- split(set[, c("SOURCE", "SPECIES.NEW", "SEQUENCE.ID", "SEQUENCE")], set$FRAGMENT)
for (s in names(sets)) {
 subs <- sets[[s]]
 subs.misu <- subs[subs$SOURCE == "MISU", ] # keep ours
 subs.misu.sp <- unique(sort(subs.misu$SPECIES.NEW))
 subs.genbank <- subs[subs$SOURCE != "MISU", ] # not ours
 subs.genbank <- subs.genbank[!subs.genbank$SPECIES %in% subs.misu.sp, ] # remove species which we already have
 subs.genbank <- subs.genbank[order(subs.genbank$SPECIES.NEW, nchar(subs.genbank$SEQUENCE), decreasing=TRUE), ] # sort by length, longer first
 subs.genbank <- subs.genbank[!duplicated(subs.genbank$SPECIES), ] # take the first
 subs <- rbind(subs.misu, subs.genbank) # join back
 subs$ID <- apply(subs[, c("SPECIES.NEW", "SEQUENCE.ID")], 1, function(.x) paste0(.x, collapse="__")) # separator between species and everything else is "__", the rest is unchanged
 subs[, c("SOURCE", "SPECIES.NEW", "SEQUENCE.ID")] <- NULL # remove columns except 2
 subs <- subs[, c("ID", "SEQUENCE")] # ID first
 Write.fasta(subs, file=paste0("20_sets/", s, ".fasta"))
}
