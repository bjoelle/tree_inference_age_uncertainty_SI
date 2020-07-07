get_tree_dist_unconstrained = function(output_folder, char) {
  
  mean_dist = 0
  
  for(i in 1:length(samp_trees)) {
    print(i)
    true_tree = samp_trees[[i]]
    attr(true_tree, "class") = c("phylo")
    
    revbayes_trees = load.trees(file.path(output_folder, paste0("output/rep_n",char,"_",i,".trees")), format = format)
    revbayes_trees = revbayes_trees[[1]]
    
    n = length(revbayes_trees)
    revbayes_trees = revbayes_trees[(round(n/10)):n]
    
    dist = c()
    for(tree in revbayes_trees) dist = c(dist, RF.dist(tree, true_tree, check.labels = TRUE, rooted=FALSE, normalize = TRUE))
    mean_dist = mean_dist + mean(dist)
  }
  
  mean_dist = mean_dist/length(samp_trees)
  return(mean_dist)
}

plot_unconstrained_priors = function(treefile, output_dir) {
  load(treefile)
  library(ape)
  library(phangorn)
  library(rwty)
  
  output_folder = paste0(output_dir, "/exp_alpha_fixed")
  exp_alpha_fixed30 = get_tree_dist_unconstrained(output_folder, char = 30) #0.8566091
  exp_alpha_fixed300 = get_tree_dist_unconstrained(output_folder, char = 300) #0.5958119
  
  output_folder = paste0(output_dir, "/exp_alpha_est")
  exp_alpha_est30 = get_tree_dist_unconstrained(output_folder, char = 30) #0.8849693
  exp_alpha_est300 = get_tree_dist_unconstrained(output_folder, char = 300) #0.6393461
  
  output_folder = paste0(output_dir, "/dirichlet")
  exp_dir30 = get_tree_dist_unconstrained(output_folder, char = 30) #0.8682971
  exp_dir300 = get_tree_dist_unconstrained(output_folder, char = 300) #0.6202797
  
  df = data.frame(char = c(30, 300, 30, 300, 30, 300), 
                  prior = c("Exp fixed rate", "Exp fixed rate", 
                            "Exp est. rate", "Exp est. rate", 
                            "Dirichlet/Gamma", "Dirichlet/Gamma"),
                  RF = c(exp_alpha_fixed30, exp_alpha_fixed300,
                         exp_alpha_est30, exp_alpha_est300, 
                         exp_dir30, exp_dir300)
  )
  
  save(df, file = "RF_values.Rdata")
  
  level_order <- c('Exp fixed rate', 'Exp est. rate', 'Dirichlet/Gamma')
  df$char = as.factor(df$char)
  
  library(ggplot2)
  
  p1 = ggplot(df, aes(x = factor(prior, level=level_order), y = RF, color = factor(char))) + geom_point(size = 3.5, shape = 17) +
    ylab("Normalized Robinson-Foulds distance") + ylim(0,1) +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), legend.position = "right") +
    scale_color_manual(values = c("gray52", "darkgoldenrod2"), name = "Character\nlength")
  
  pdf("tree_prior.pdf", width = 7, height = 6)
  p1
  dev.off()
}