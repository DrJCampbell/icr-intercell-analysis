# icr-intercell-analysis

A collection of R functions used to run analyses of Intercell data sets. Typical usage would be:

```
setwd("/Users/jamesc/path/to/analysis/dir/")
source("/Users/jamesc/path/to/this/git/repo/0.library_functions.R")

# Define the input files
rnai_file <- "/Users/jamesc/path/to/rnai/data/zscores.txt"
tissues_file <- "tissues.txt"
func_muts_file <- "functional_mutations.txt"
all_muts_file <- "all_mutations.txt"
mut_classes_file <- "mutation_classes_150225.txt"

rnai_combmuts <- read_rnai_mutations(
	rnai_file=rnai_file,
	func_muts_file=func_muts_file,
	all_muts_file=all_muts_file,
	mut_classes_file=mut_classes_file,
	tissues_file=tissues_file
	)

run_intercell(
	x=rnai_combmuts,
	qnorm=TRUE,
	analysis_id="rnai_combmuts_analysis"
	)
```
