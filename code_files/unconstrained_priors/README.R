#This folder contains the files used to explore the impact of the branch length prior in unconstrained (i.e. non-time calibrated) tree inference.

# To create the Rev scripts used in the inference
# bash rev_files.sh

# To summarize results and make the plot
source("priors_figure.R")
treefile = "path/to/tree/file" #highfs or lowfs, trees are the same for all charlengths
output_folder = "path/to/output/files"
plot_unconstrained_priors(treefile, output_folder)