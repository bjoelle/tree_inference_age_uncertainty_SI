source("crinoids_log_analysis.R")
source("crinoids_tree_analysis.R")
source("read.incomplete.nexus.R")
source("run_analysis.R")

# Run the full analysis of BEAST2/RevBayes output on a simulated dataset
# output will be considered BEAST2 unless the folder name contains "unconstrained"
treefiles = "path/to/folder/of/datasets/stored/as/RData"
outfolder = "path/to/folders/of/BEAST2/logfiles"
resultfolder = "path/to/store/results"
prior = TRUE # true if folder paths contain the prior name, false otherwise
run_analysis(outfolder, treefiles, resultfolder, prior = prior)

# Some runs were skipped in the final analysis because they did not converge
# A list of these can be found in "skip.RData", where skip_300 and skip_3000 correspond
# respectively to the highfs_charlength300 and highfs_charlength3000 runs
load("skip.RData")