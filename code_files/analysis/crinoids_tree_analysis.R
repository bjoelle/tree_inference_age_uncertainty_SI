FBD_allnodes_divtimes = function(treefile, outfolder, tmp_file = "tmp_analysis.RData") {
  library(ape)
  library(coda)
  load(treefile)
  
  names = c("correct_ages", "median_ages", "interval_ages", "random_ages", "norm_intvl_ages")
  if(file.exists(tmp_file)) {
    load(tmp_file)
    analysis = get("analysis")
    initi = length(analysis[[1]])+1
  } else {
    analysis = list()
    initi = 1
  }
  if(initi > length(trees)) return(analysis)
  
  for(i in initi:length(trees)) {
    print(i)
    mrca_list = list(ages = c(), tips = list())
    used_nodes = c()
    ages = n.ages(samp_trees[[i]])
    ages = ages + offsets[[i]]$correct_ages
    used_tips = length(samp_trees[[i]]$tip.label)
    
    all_mrca = mrca(samp_trees[[i]])
    for(t1 in 1:(used_tips-1)) {
      idx1 = which(rownames(all_mrca) == samp_trees[[i]]$tip.label[t1])
      for(t2 in (t1+1):used_tips) {
        idx2 = which(colnames(all_mrca) == samp_trees[[i]]$tip.label[t2])
        mrca = all_mrca[idx1,idx2]
        if(!mrca %in% used_nodes) {
          used_nodes = c(used_nodes, mrca)
          mrca_list$tips[[length(mrca_list$tips)+1]] = c(samp_trees[[i]]$tip.label[t1], samp_trees[[i]]$tip.label[t2])
          mrca_list$ages = c(mrca_list$ages, ages[mrca])
        }
      }
    }
    
    for (nm in names) {
      print(nm)
      beast_trees = read.incomplete.nexus(file.path(outfolder, paste0("FBD_", nm,"_",i,".trees")))
      beast_trees = lapply(beast_trees, function(x) {x$tip.label = attr(beast_trees,"TipLabel"); x})
      nt = length(beast_trees)
      
      if(nm == "interval_ages" || nm == "norm_intvl_ages") {
        offs = read.table(file.path(outfolder, paste0("FBD_", nm,"_",i,".log")),header = T,comment.char = "#")$offset
        no = length(offs)
        if(nt == no*2 || nt == no*2-1) {
          for(iii in length(beast_trees):1) if(iii %% 2 == 0) beast_trees = beast_trees[-iii]
          nt = length(beast_trees)
        }
        #if(nt != no) browser()
        offs = offs[(round(no/10)+1):no]
      } else offs = rep(offsets[[i]][[nm]], length(beast_trees))
      
      nt = length(beast_trees)
      beast_trees = beast_trees[(round(nt/10)+1):nt]
      
      tmp = list()
      for(mr in mrca_list$tips) tmp[[paste0(mr[1], ".", mr[2])]] = rep(0, length(beast_trees))
      for (bti in 1:length(beast_trees)) {
        ages = n.ages(beast_trees[[bti]])
        for(mr in mrca_list$tips) {
          t = ages[getMRCA(beast_trees[[bti]], mr)] + offs[bti]
          tmp[[paste0(mr[1], ".", mr[2])]][bti] = t
        }
      }
      rel_error = c()
      cover = wdth = c()
      for(j in 1:length(mrca_list$ages)) {
        rel_error = c(rel_error, abs(median(tmp[[paste0(mrca_list$tips[[j]][1], ".", mrca_list$tips[[j]][2])]]) 
                                     - mrca_list$ages[j])/mrca_list$ages[j])
        if(any(is.na(rel_error))) browser()
        HPD = HPDinterval(as.mcmc(tmp[[paste0(mrca_list$tips[[j]][1], ".", mrca_list$tips[[j]][2])]]))
        cover = c(cover, (mrca_list$ages[j] >= HPD[1] - 1e-6) && (mrca_list$ages[j] <= HPD[2] + 1e-6))
        wdth = c(wdth, (HPD[2] - HPD[1])/mrca_list$ages[j])
      }
      
      analysis[[paste0(nm,".extmrca.divtime.relative_error")]] = c(analysis[[paste0(nm,".extmrca.divtime.relative_error")]], mean(rel_error))
      analysis[[paste0(nm,".extmrca.divtime.coverage")]] = c(analysis[[paste0(nm,".extmrca.divtime.coverage")]], sum(cover)/length(cover))
      analysis[[paste0(nm,".extmrca.divtime.HPD_width")]] = c(analysis[[paste0(nm,".extmrca.divtime.HPD_width")]], mean(wdth))
    }
    save(analysis,file=tmp_file)
  }
  analysis
}

n.ages = function(tree){
  
  node.ages <- TreeSim::getx(tree, sersampling = 1)[1:(tree$Nnode+length(tree$tip))]
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))
  
  # adding possible offset if tree fully extinct
  if(!is.null(tree$root.time)) node.ages = node.ages + tree$root.time - max(node.ages)
  
  return(node.ages)
}

get_RF_dist = function(treefile, output_folder, save_folder) {
  library(ape)
  library(phangorn)
  
  names = c("correct_ages", "median_ages", "interval_ages", "random_ages", "norm_intvl_ages")
  load(treefile)
  mean_dist = list()
  
  for(i in 1:length(samp_trees)) {
    print(i)
    true_tree = samp_trees[[i]]
    attr(true_tree, "class") = c("phylo")
    
    for(name in names) {
      print(name)
      beast_trees = read.incomplete.nexus(file.path(output_folder, paste0("FBD_", name,"_",i,".trees")))
      beast_trees = lapply(beast_trees, function(x) {x$tip.label = attr(beast_trees,"TipLabel"); x})
      n = length(beast_trees)
      beast_trees = beast_trees[(round(n/10)):n]
      
      dist = c()
      for(tree in beast_trees) dist = c(dist, RF.dist(tree, true_tree, check.labels = TRUE, rooted=FALSE, normalize = TRUE))
      mean_dist[[name]] = c(mean_dist[[name]], mean(dist))
    }
  }
  
  save.file = paste0(save_folder, basename(output_folder), "_mean_norm_dist.RData")
  assign(paste0(basename(output_folder),"_mean_norm_dist"), mean_dist)
  save(list = c(paste0(basename(output_folder), "_mean_norm_dist")), file = save.file)
}

get_RF_dist_unconstrained = function(treefile, output_folder, save_folder) {
  library(ape)
  library(phangorn)
  library(rwty)
  
  load(treefile)
  mean_dist = c()
  
  for(i in 1:length(samp_trees)) {
    print(i)
    true_tree = samp_trees[[i]]
    attr(true_tree, "class") = c("phylo")
    
    revbayes_trees = load.trees(file.path(output_folder, paste0("output/rep_",i,".trees")), format = "revbayes")
    revbayes_trees = revbayes_trees[[1]]
    n = length(revbayes_trees)
    revbayes_trees = revbayes_trees[(round(n/10)):n]
    
    dist = c()
    for(tree in revbayes_trees) dist = c(dist, RF.dist(tree, true_tree, check.labels = TRUE, rooted=FALSE, normalize = TRUE))
    mean_dist = c(mean_dist, mean(dist))
  }
  
  save.file = paste0(save_folder, basename(output_folder), "_mean_norm_dist.RData")
  assign(paste0(basename(output_folder),"_mean_norm_dist"), mean_dist)
  save(list = c(paste0(basename(output_folder), "_mean_norm_dist")), 
       file = save.file)
}
