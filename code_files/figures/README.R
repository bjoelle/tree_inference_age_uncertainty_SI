source("figures.R")
source("fig_empirical.R")

# To make a plot of all prior comparison results
outfolder = "path/to/analysis/files/as/RData"
plotfolder = "path/to/store/figures"
make_figures_priors(outfolder, plotfolder)

# To make a plot of all fossil and character sampling comparison results
outfolder = "path/to/analysis/files/as/RData"
plotfolder = "path/to/store/figures"
make_figures_sampling(outfolder, plotfolder)

# To make a boxplot comparing estimates of FBD parameters on brachiopods
outfolder = "path/to/folder/with/BEAST2/logfiles"
plotfile = "path/to/output/figure"
brachiopods_boxplot(outfolder, plotfile)