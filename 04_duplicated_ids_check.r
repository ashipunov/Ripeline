## this script checks if "dna" table contains duplicated ids
## it is sometimes possible because some IDs are compound and amount of data large
## this script is not a part of main pipeline, it must be run separately when data changes

dna <- read.table("_kubricks_dna.txt", h=TRUE, sep="\t", as.is=TRUE)
dna <- dna[dna$SELECT == 1, ]

ids <- strsplit(dna$SEQUENCE.ID, ", *")
for (i in 1:length(ids)) ids[[i]] <- paste(ids[[i]], dna$FRAGMENT[i], sep="_")
ids <- unlist(ids)
ids[duplicated(ids)]
