run_analysis = function(outfolder, datafolder, savefolder, prior = T, add = "_seed459", overwrite = F) {
  dirs = list.dirs(outfolder, recursive = F)
  for (d in dirs) {
    treef = strsplit(basename(d),"_")[[1]]
    if(prior) treef = treef[-length(treef)]
    treef = paste0(treef, collapse = "_")
    treef = paste0(treef, add, ".RData")
    treef = file.path(datafolder, treef)
    if(!grepl("unconstrained", d)) {
      if(overwrite || !file.exists(file.path(savefolder, paste0(basename(d),"_summary.RData")))) 
         FBD_analysis(treef, d, save_folder = savefolder)
      if(overwrite || !file.exists(file.path(savefolder, paste0(basename(d),"_mean_norm_dist.RData")))) 
        get_RF_dist(treef, d, savefolder)
    }
    else {
      if(overwrite || !file.exists(file.path(savefolder, paste0(basename(d),"_mean_norm_dist.RData"))))
        get_RF_dist_unconstrained(treef, d, savefolder)
    }
  }
}