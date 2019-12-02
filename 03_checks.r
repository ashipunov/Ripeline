## this script checks if "dna" and "treesp" tables are concerted with the main "sp" table
## some minor (but explainable) differences are typically allowed
## this script is not a part of main pipeline, it must be run separately when data changes

## read input (file names start with "_") tables
dna <- read.table("_kubricks_dna.txt", sep="\t", h=TRUE, as.is=TRUE)
treesp <- read.table("_kubricks_treesp.txt", sep="\t", h=TRUE, as.is=TRUE)
sp <- read.table("_kubricks_sp.txt", sep="\t", h=TRUE, as.is=TRUE)

## check against "sp" table
cat("kubricks dna not in sp:\n", paste0(sort(setdiff(dna$SPECIES.NEW, sp$SPECIES)), "\n"), "\n")
cat("kubricks treesp not in sp:\n", paste0(sort(setdiff(treesp$SPECIES, sp$SPECIES)), "\n"), "\n")
