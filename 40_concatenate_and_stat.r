## this script takes alignments, concatenates sequences (i.e., makes super-matrices) and outputs some statistics

library(shipunov)

## we will need DNA database again to access all species and IDs
sequences <- read.table("_kubricks_dna.txt", sep="\t", h=TRUE, as.is=TRUE)
## we will need "treesp" to determine outgroups
treesp <- read.table("_kubricks_treesp.txt", sep="\t", h=TRUE, as.is=TRUE)
## for subsequent analyses, outgroups are better to place in the beginning of file so we need to know what are our outgroups
outgroups <- treesp$SPECIES.NEW[treesp$TYPE == "outgroup"]

## now load the actual DNA data
## define fragments
fragments <- c("abcd", "efgh") # fragile (this means that if you must control this when e.g. fragments change)
## load alignments per fragment
dna <- structure(vector("list", length(fragments)), names=fragments) # fragile
for (f in fragments) {
 dna[[f]] <- Read.fasta(paste0("32_alignments_trimmed_gapcoded/", f, ".fasta"))
 sp.id <- do.call(rbind, strsplit(dna[[f]]$name, "__")) # we use here that species name reparated with everything else by _double_ underscore
 dna[[f]] <- data.frame(SP=sp.id[, 1], ID=sp.id[, 2], SEQ=dna[[f]][, 2], stringsAsFactors=FALSE)
 dna[[f]] <- dna[[f]][order(nchar(gsub("-", "", dna[[f]]$SEQ)), decreasing=TRUE), ] # longest first---it is useful for concatenation
 dna[[f]] <- dna[[f]][order(dna[[f]]$SP), ] # by species, whithin species longest still first because order() does not rearrange ties
}

## now make the list of all IDs with their taxonomic names
## the main problem now is that some IDs are _compound_ and must be split
sp.ids <- do.call(rbind, c(lapply(dna, "[", 1:2), make.row.names=FALSE))
tmp1 <- suppressWarnings(do.call(rbind, strsplit(sp.ids$ID, ", "))) # as some IDs are compound, and some not, strsplit() will issue warning and it is better to suppress it
tmp2 <- unique(sort(paste(c(tmp1), sp.ids$SP, sep="__"))) # for simpler sort / unique, IDs and names join with "__" again
ids.sp <- do.call(rbind, strsplit(tmp2, "__")) # then we split them back again. At the end, "ids.sp" will contain the unique IDs and their (species) names

## ===
## ===
## ===

## STEP 1. STRICT concatenation: fragments merge on "our" ids. In other words, only our (locally made) sequences with the same ID will be concatenated.

## determine which IDs are "ours"
tmp1 <- sequences[sequences$SOURCE == "MISU", "SEQUENCE.ID"]
tmp2 <- suppressWarnings(do.call(rbind, strsplit(tmp1, ", ")))
ourids.all <- unique(sort(tmp2))
ourids <- ids.sp[ids.sp[, 1] %in% ourids.all, ]

## make the dataframe to keep concatenated
strict <- data.frame(ourids, stringsAsFactors=FALSE)
names(strict) <- c("ID", "SP")

## reorder to place outgroups first, other sorted
tmp1 <- strict[strict$SP %in% outgroups, ]
tmp2 <- strict[!strict$SP %in% outgroups, ]
tmp2 <- tmp2[order(tmp2$SP), ]
strict <- rbind(tmp1, tmp2)

## concatenate!
for (n in fragments) { # per each fragment
 strict[[n]] <- NA # empty column, filled with NAs
 for (i in 1:nrow(strict)) {
  ii <- dna[[n]]$SEQ[grepl(strict$ID[i], dna[[n]]$ID)] # search for sequence which ID corresponds with ID in the current row of "strict"
  strict[[n]][i] <- ifelse(length(ii) == 0, NA, ii) # if by any chance "ii" has length > 1, returns the first element; if "ii" is empty, make NA
 }
}
## at the end, "strict" has columns for ID, SP(ecies) and one column for each fragment. To finalize contatenation, we need to join fragment columns (see below)

## ===
## ===
## ===

## STEP 2. SEMISTRICT concatenation: first, add new species which were not in "our" ids; then fill gaps with sequences from the same species and keep these IDs

## make a copy
semistrict <- strict

## add species which are not "ours"
ssp <- sort(unique(semistrict$SP))
allsp <- sort(unique(sp.ids$SP))
addsp <- allsp[!allsp %in% ssp]
if (length(addsp) > 0) { # precaution for the case when all species are locally sequenced
semistrict[(nrow(semistrict) + 1) : (nrow(semistrict) + length(addsp)), ] <- NA # fill new rows with NAs
semistrict[nrow(semistrict) : (nrow(semistrict) - (length(addsp) - 1)), "SP"] <- addsp # number of rows changed, hence new calculation; fill new part of SP column with "addsp"
}

## reorder again to place outgroups first, other sorted
tmp1 <- semistrict[semistrict$SP %in% outgroups, ]
tmp2 <- semistrict[!semistrict$SP %in% outgroups, ]
tmp2 <- tmp2[order(tmp2$SP), ]
semistrict <- rbind(tmp1, tmp2)

## make new columns to keep IDs of will-be-added sequences
names.tmp3 <- c(paste0(fragments, ".ID"))
tmp3 <- structure(data.frame(matrix(NA, nrow=nrow(semistrict), ncol=length(names.tmp3))), names=names.tmp3)
semistrict <- cbind(semistrict, tmp3)
semistrict <- semistrict[, c("SP", "ID", paste0(fragments, ".ID"), fragments)]

## concatenate!
for(n in fragments) { # per fragment
 d <- dna[[n]]
 f.id <- paste0(n, ".ID") # make the fragment ID to determine where to place it
 for (i in 1:nrow(semistrict)) {
  if (is.na(semistrict[i, n])) { # if there is a gap
  dd <- d$SEQ[grepl(semistrict$SP[i], d$SP)] # search for sequence of this fragment with same species name
  ii <- d$ID[grepl(semistrict$SP[i], d$SP)] # analogously, search for ID
  if (length(dd) > 0 & length(ii) > 0) { # if anything found
   semistrict[i, n] <- dd[1] # then select the first (i.e., longest because they were sorted this way) element
   semistrict[i, f.id] <- ii[1] # and select its ID
   } else { # if nothing found
   semistrict[i, f.id] <- "-" # put there the dash
   }
 }}
}
## at the end, "semistrict" has columns for our SP, our ID, then columns for added IDs (one per fragment), then columns for each fragment. If, for example, there are 3 fragments, "semistrict" will have 8 columns.

## ===
## ===
## ===

## STEP 3. Output statistics and write files
## we will output only semistrict concatenation

## create new IDs which will hshow the origin of tree terminal
semistrict.f.ids <- c(paste0(fragments, ".ID"))
ids.tmp0 <- semistrict[, c("SP", "ID", semistrict.f.ids)]
ids.tmp1 <- apply(ids.tmp0[, c("SP", "ID")], 1, function(.x) paste(.x, collapse="__")) # as usual, species separated from ID with "__"
for (i in semistrict.f.ids) ids.tmp0[, i] <- ifelse(is.na(ids.tmp0[, i]), "", paste(sub(".ID", "", i), ids.tmp0[, i])) # add name of fragment to sequence ID
ids.tmp2 <- apply(ids.tmp0[, semistrict.f.ids], 1, function(.x) paste(.x, collapse=" ")) # join them all together
ids.tmp2 <- gsub("(^ +)|( +$)", "", ids.tmp2) # then clean ...
ids.tmp2 <- gsub("( +)", " ", ids.tmp2)
ids.tmp <- paste(ids.tmp1, ids.tmp2)
ids.tmp <- gsub(" +", "_", ids.tmp)
ids.tmp <- gsub(",", "", ids.tmp)
ids.tmp <- gsub("_$", "", ids.tmp)
ids.tmp <- gsub("_NA_", "_", ids.tmp) # ... cleaning done
## at the end, new ID will start with species name, then "our" ID (if any), then name of fragment for added "alien" sequence, then its ID or "-" if fragment absent (please look on trees to understand this system better)

## remove duplicated sequences, they appear partly because some sequences in DNA database had compound IDs
semistrict.f <- fragments
groups <- as.numeric(as.factor(apply(semistrict[, c("SP", semistrict.f)], 1, paste, collapse="")))
## since some sequences were essentially removed (even if they are full duplicates), it is good to kee their IDs
comb.ids <-  aggregate(ids.tmp, by=list(groups), paste, collapse=" + ")
semistrict <- semistrict[!duplicated(groups), ] # remove duplicates in main table
ids <- ids.tmp[!duplicated(groups)] # remove duplicates in joint ids
comb.ids <- comb.ids[groups[!duplicated(groups)], 2] # the same order as in deduplicated data
## this table contains combined IDs and is useful to understand which IDs are omitted in resulted trees
write.table(data.frame(ids, comb.ids), file="40_concatenated/semistrict_comb_ids.txt", sep="\t", quote=FALSE, row.names=FALSE)

## output statistics:
semistrict.m <- data.frame(semistrict[, 1:2], apply(semistrict[, semistrict.f], 2, function(.y) ifelse(is.na(.y), 0, 1)))
semistrict.m$COUNT <- rowSums(semistrict.m[, semistrict.f])
## what missed after concatenation
semistrict.m[order(semistrict.m$COUNT, decreasing=TRUE), ]
## how many ids in each missing state
table(apply(semistrict[, semistrict.f], 1, function(.y) sum(!is.na(.y))))
## how many ids missed in each fragment
apply(semistrict[, semistrict.f], 2, function(.x) paste0(sum(is.na(.x)), "/", length(.x)))
## how many sequences our and total
all.ids <- unlist(semistrict[, c("ID", paste0(semistrict.f, ".ID"))])
## our
sum(all.ids %in% ourids.all & !is.na(all.ids) & all.ids != "-")
## total
sum(!is.na(all.ids) & all.ids != "-")
## number of potentially informative characters and total number of characters per fragment
for (f in semistrict.f) {
cat(f)
tmp <- semistrict[, f]
tmp <- tmp[!is.na(tmp)]
tmp.table <- do.call(rbind, strsplit(tmp, split=""))
tmp.table[tmp.table == "-" | tmp.table == "N"] <- NA
cat("\tall chars ", ncol(tmp.table))
tmp.itic <- apply(tmp.table, 2, Is.tax.inform.char)
cat("\tpotentially taxonomically informative chars ", sum(tmp.itic), "\n")
}

## to make FASTA, we need to replace NA with long strings filled with "N"
for (i in semistrict.f) semistrict[, i] <- ifelse(is.na(semistrict[, i]), paste(rep("N", max(nchar(semistrict[, i]), na.rm=TRUE)), collapse=""), semistrict[, i])
## now join all sequences columns into one
seq <- apply(semistrict[, semistrict.f], 1, function(.x) paste(.x, collapse=""))
semistrict.new <- data.frame(ID=ids, SEQ=seq, stringsAsFactors=FALSE)
Write.fasta(semistrict.new, file=paste0("40_concatenated/", "semistrict.fasta"))
