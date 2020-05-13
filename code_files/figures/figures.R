make_figures_priors = function(outfolder, plotfolder) {
  plot_tree_distances_priors(outfolder, T, plotfile = file.path(plotfolder, "highfs_norm_tree_distance.pdf"))
  plot_tree_distances_priors(outfolder, F, plotfile = file.path(plotfolder, "lowfs_norm_tree_distance.pdf"))
  plot_comparison_priors(outfolder, plotfolder = plotfolder)
}

make_figures_sampling = function(outfolder, plotfolder) {
  plot_tree_distances_sampling(outfolder, T, plotfile = file.path(plotfolder, "highfs_norm_tree_distance.pdf"))
  plot_tree_distances_sampling(outfolder, F, plotfile = file.path(plotfolder, "lowfs_norm_tree_distance.pdf"))
  plot_comparison_sampling(outfolder, plotfolder = plotfolder)
}

plot_tree_distances_priors = function(outfolder, highfs, priors = c("nopriors", "clockprior","correctclockprior"),
                                      plotfile = NULL) {
  library(ggplot2)
  ages = c("*Interval ages", "*Symmetric ages", "Correct age", "Median age", "Random age")
  cats = c("interval_ages", "norm_intvl_ages", "correct_ages", "median_ages", "random_ages")
  df = NULL
  
  for(prior in priors) {
    file_name = paste0(outfolder, "/", 
                       if(highfs) "highfs_charlength300_" else "lowfs_charlength30_", prior, "_mean_norm_dist.RData")
    load(file_name)
    tree_dist = get(tools::file_path_sans_ext(basename(file_name)))
    for(i in 1:length(cats))  df = rbind(df, data.frame(Age = ages[i], mean_RF_distance = mean(tree_dist[[cats[i]]], na.rm = T),
                                                        sd_RF_distance = sd(tree_dist[[cats[i]]], na.rm = T), Prior = prior))
  }
  
  unconstrained_file_name = paste0(outfolder, "/", 
                                   if(highfs) "highfs_charlength300_" else "lowfs_charlength30_", "unconstrained_mean_norm_dist.RData")
  load(unconstrained_file_name)
  tree_dist = get(tools::file_path_sans_ext(basename(unconstrained_file_name)))
  tmp_df = data.frame(Age = "No ages (unconstrained)", mean_RF_distance = mean(tree_dist, na.rm = T), 
                      sd_RF_distance = sd(tree_dist, na.rm = T), Prior = "unconstrained")
  df = rbind(df, tmp_df)
  
  cbPalette <- c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  if(highfs) s = 16 else s = 17
  df$Age = factor(df$Age, levels = c(ages, "No ages (unconstrained)"))
  
  df$lower = sapply(df$mean_RF_distance-df$sd_RF_distance, function(t) max(0,t))
  
  p2 = ggplot(df, mapping=aes(x=Age, y=mean_RF_distance, colour=factor(Prior))) + xlab(NULL) + 
    geom_pointrange(aes(ymin=lower, ymax=mean_RF_distance+sd_RF_distance), size = 1, shape = s, position = position_dodge(width=0.75)) +
    xlim("*Interval ages", "*Symmetric ages", "Correct age", "Median age", "Random age", "No ages (unconstrained)") +
    ylab("Normalized Robinson-Foulds distance") + coord_cartesian(ylim = c(0,1)) + #geom_jitter(size = 3.5, shape = s, width=0.25, height=0) +
    theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), legend.position = "bottom",
          legend.direction = "vertical") + scale_color_manual(values=cbPalette, name = "Clock rate prior", breaks = c(priors, "unconstrained"),
                                                              labels = c("Unbounded uniform prior", "Lognormal prior with rate mean 0.2/Myr", "Lognormal prior with rate mean 0.033/Myr", "Unconstrained (no priors)"))
  
  if(is.null(plotfile)) show(p2)
  else ggsave(plotfile,width = 5,height = 7)
}

plot_tree_distances_sampling = function(outfolder, highfs, charlengths = c("30", "300", "3000"),
                                        plotfile = NULL) {
  library(ggplot2)
  ages = c("*Interval ages", "*Symmetric ages", "Correct age", "Median age", "Random age",  "No ages (unconstrained)")
  cats = c("interval_ages", "norm_intvl_ages", "correct_ages","median_ages", "random_ages", "unconstrained")
  tree_dist = list()
  df = NULL
  
  for(cl in charlengths) {
    file_name = paste0(outfolder, "/", 
                       if(highfs) "highfs_" else "lowfs_", "charlength", cl, "_mean_norm_dist.RData")
    load(file_name)
    
    unconstrained_file_name = paste0(outfolder, "/",
                                     if(highfs) "highfs_" else "lowfs_", "charlength", cl, 
                                     "_unconstrained_mean_norm_dist.RData")
    load(unconstrained_file_name)
    
    tree_dist[[cl]] = c(get(tools::file_path_sans_ext(basename(file_name))), list(get(tools::file_path_sans_ext(basename(unconstrained_file_name)))))
    names(tree_dist[[cl]])[length(cats)] = "unconstrained"
    for(i in 1:length(cats)) {
      df = rbind(df, data.frame(Age = ages[i], mean_RF_distance = mean(tree_dist[[cl]][[cats[i]]], na.rm = T), 
                                sd_RF_distance = sd(tree_dist[[cl]][[cats[i]]], na.rm = T), Char_length = cl))
    }
  }
  
  df$Age = factor(df$Age, levels = ages)
  cbPalette <- c("#999999", "#E69F00", "#56B4E9")
  if(highfs) s = 16 else s = 17
  
  df$lower = sapply(df$mean_RF_distance-df$sd_RF_distance, function(t) max(0,t))
  
  p2 = ggplot(df, mapping=aes(x=Age, y=mean_RF_distance, colour=factor(Char_length))) + xlab(NULL) + 
    geom_pointrange(aes(ymin=lower, ymax=mean_RF_distance+sd_RF_distance), size = 1, shape = s, position = position_jitter(width=0.25, height=0)) +
    xlim("*Interval ages", "*Symmetric ages", "Correct age", "Median age", "Random age", "No ages (unconstrained)") +
    ylab("Normalized Robinson-Foulds distance") + coord_cartesian(ylim = c(0,1)) + #geom_jitter(size = 3.5, shape = s, width=0.25, height=0) +
    theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), legend.position = "bottom",
          legend.direction = "vertical") + scale_color_manual(values=cbPalette, name = "Morphological character length")
  if(is.null(plotfile)) show(p2)
  else ggsave(plotfile, width = 5, height = 7)
}

plot_comparison_priors = function(outfolder, priors = c("nopriors", "clockprior","correctclockprior"),
                                  names = c("divtime", "diversificationRateFBD","turnoverFBD", "clockRate"),
                                  plotfolder = NULL) {
  library(ggplot2)
  library(gridExtra)
  plotnames = c("Divergence time", "Diversification rate", "Turnover", "Clock rate")
  names(plotnames) = c("extmrca.divtime", "diversificationRateFBD","turnoverFBD", "clockRate")
  ages = c("Correct age", "Median age", "*Interval ages", "Random age", "*Symmetric ages")
  cats = c("correct_ages", "median_ages", "interval_ages", "random_ages", "norm_intvl_ages")
  
  for(name in names) {
    if(name == "divtime") name = "extmrca.divtime"
    assign(paste0(name, "_rel_error"), NULL)
    assign(paste0(name, "_coverage"), NULL)
  }
  
  for(prior in priors) {
    for(type in c("highfs_charlength300", "lowfs_charlength30")) {
      file = paste0(outfolder, "/", type, "_", prior, "_analysis.RData")
      
      load(file)
      fullfile = get(tools::file_path_sans_ext(basename(file)))
      for(name in names) {
        full = fullfile[[name]]
        if(name == "divtime") name = "extmrca.divtime"
        
        tmp_df = NULL
        for(i in 1:length(cats)) {
          vect = full[[paste0(cats[i],".",name,".relative_error")]]
          if(sum(is.infinite(vect)) != 0) vect = vect[-which(is.infinite(vect))]
          tmp_df = rbind(tmp_df, data.frame(Age = ages[i], mean_Relative_error = mean(vect, na.rm = T),
                                            sd_Relative_error = sd(vect, na.rm = T), Prior = prior, Type = type))
        }
        assign(paste0(name, "_rel_error"), rbind(get(paste0(name, "_rel_error")), tmp_df))
        
        tmp_df = data.frame(Age = ages, Coverage = unlist(lapply(cats, function(x) {mean(full[[paste0(x,".",name,".coverage")]], na.rm = T)})),
                            Prior = rep(prior, length(cats)), Type = rep(type, length(cats)))
        assign(paste0(name, "_coverage"), rbind(get(paste0(name, "_coverage")), tmp_df))
      }
    }
  }
  
  cbPalette <- c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  for(name in names) {
    if(name == "divtime") name = "extmrca.divtime"
    
    tmp_df = get(paste0(name, "_rel_error"))
    tmp_df$lower = sapply(tmp_df$mean_Relative_error-tmp_df$sd_Relative_error, function(t) max(0,t))
    
    p_rel_error = ggplot(tmp_df, mapping = aes(x=Age, y=mean_Relative_error, colour=factor(Prior), shape = factor(Type))) +
      geom_pointrange(aes(ymin=lower, ymax=mean_Relative_error+sd_Relative_error), 
                      size = 0.75, position = position_dodge(width=0.75)) +
      xlim("*Interval ages", "*Symmetric ages", "Correct age", "Median age", "Random age") +
      coord_cartesian(ylim = c(0,1)) + xlab(plotnames[name]) + ylab("Relative error") #+ geom_jitter(aes(shape=factor(Type)), size = 4, width = 0.3, height = 0)
    
    if(name == "turnoverFBD" | name == "clockRate") {
      p_rel_error = p_rel_error +
        theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), 
              legend.position = "bottom", legend.direction = "vertical", legend.box = "vertical") +
        scale_color_manual(values=cbPalette, name = "Clock rate prior", breaks = priors, 
                           labels = c("Unbounded uniform prior", "Lognormal prior (median not equal to true rate)", "Lognormal prior (median equal to true rate)")) + 
        scale_shape_manual(values = c(16, 17), name = "Fossil sampling and character length", breaks = c("highfs_charlength300", "lowfs_charlength30"), 
                           labels = c("High sampling and character length 300", "Low sampling and character length 30"))
      
      if(is.null(plotfolder)) show(p_rel_error)
      else ggsave(paste0(plotfolder, "/", if(name == "extmrca.divtime") "divtime" else name, "_rel_error.pdf"),width = 5,height = 7.25)
    } else {
      p_rel_error = p_rel_error +
        theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), legend.position = "none") +
        scale_color_manual(values=cbPalette) + scale_shape_manual(values = c(16, 17))
      
      if(is.null(plotfolder)) show(p_rel_error)
      else ggsave(paste0(plotfolder, "/", if(name == "extmrca.divtime") "divtime" else name, "_rel_error.pdf"),width = 5,height = 5)
    }
    
    p_coverage = ggplot(get(paste0(name, "_coverage")), mapping = aes(x=Age, y=Coverage, colour=factor(Prior), shape = factor(Type))) +
      xlim("*Interval ages", "*Symmetric ages", "Correct age", "Median age", "Random age") +
      ylim(0,1) + xlab(plotnames[name]) + ylab("Coverage") + geom_point(size = 3, position = position_dodge(width=0.75)) +
      theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), legend.position = "none") +
      scale_color_manual(values=cbPalette) + scale_shape_manual(values = c(16, 17))
    
    if(is.null(plotfolder)) show(p_coverage)
    else ggsave(paste0(plotfolder, "/", if(name == "extmrca.divtime") "divtime" else name, "_coverage.pdf"),width = 5,height = 5)
  }
}

plot_comparison_sampling = function(outfolder, charlengths = c("30", "300", "3000"),
                                    names = c("divtime", "diversificationRateFBD","turnoverFBD", "clockRate"),
                                    plotfolder = NULL) {
  library(ggplot2)
  plotnames = c("Divergence time", "Diversification rate","Turnover", "Clock rate")
  names(plotnames) = c("extmrca.divtime", "diversificationRateFBD","turnoverFBD", "clockRate")
  
  cats = c("correct_ages", "median_ages", "interval_ages", "random_ages", "norm_intvl_ages")
  ages = c("Correct age", "Median age", "*Interval ages", "Random age", "*Symmetric ages")
  
  for(name in names) {
    if(name == "divtime") name = "extmrca.divtime"
    assign(paste0(name, "_rel_error"), NULL)
    assign(paste0(name, "_coverage"), NULL)
  }
  for(cl in charlengths) {
    for(type in c("highfs", "lowfs")) {
      file = paste0(outfolder, "/", type, "_charlength", cl, "_analysis.RData")
      load(file)
      fullfile = get(tools::file_path_sans_ext(basename(file)))
      
      for(name in names) {
        full = fullfile[[name]]
        if(name == "divtime") name = "extmrca.divtime"
        
        tmp_df = NULL
        for(i in 1:length(cats)) {
          vect = full[[paste0(cats[i],".",name,".relative_error")]]
          tmp_df = rbind(tmp_df, data.frame(Age = ages[i], mean_Relative_error = mean(vect, na.rm = T),
                                            sd_Relative_error = sd(vect, na.rm = T), Char_length = cl, Type = type))
        }
        assign(paste0(name, "_rel_error"), rbind(get(paste0(name, "_rel_error")), tmp_df))
        
        tmp_df = data.frame(Age = ages, Coverage = unlist(lapply(cats, function(x) {mean(full[[paste0(x,".",name,".coverage")]], na.rm = T)})),
                            Char_length = rep(cl, length(cats)), Type = rep(type, length(cats)))
        assign(paste0(name, "_coverage"), rbind(get(paste0(name, "_coverage")), tmp_df))
      }
    }
  }
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  for(name in names) {
    if(name == "divtime") name = "extmrca.divtime"
    
    tmp_df = get(paste0(name, "_rel_error"))
    tmp_df$lower = sapply(tmp_df$mean_Relative_error-tmp_df$sd_Relative_error, function(t) max(0,t))
                 
    p_rel_error <- ggplot(tmp_df, mapping = aes(x=Age, y=mean_Relative_error, colour=factor(Char_length), shape = factor(Type))) +
      geom_pointrange(aes(ymin = lower, ymax=mean_Relative_error+sd_Relative_error), 
                      size = 0.75, position = position_dodge(width=0.75)) +
      xlim("*Interval ages", "*Symmetric ages", "Correct age", "Median age", "Random age") +
      coord_cartesian(ylim = c(0,1)) + xlab(plotnames[name]) + ylab("Relative error") #+ geom_jitter(aes(shape=factor(Type)), size = 4, width=0.3, height=0)
    
    if(name == "turnoverFBD" | name == "clockRate") {
      p_rel_error = p_rel_error + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), 
                                        legend.position = "bottom", legend.direction = "vertical", legend.box = "vertical") +
        scale_color_manual(values = cbPalette, name = "Morphological character length", breaks = charlengths, labels = charlengths) +
        scale_shape_manual(values = c(16, 17), name = "Fossil sampling", breaks = c("highfs", "lowfs"), 
                           labels = c("High sampling", "Low sampling"))
      if(is.null(plotfolder)) show(p_rel_error)
      else ggsave(paste0(plotfolder, "/", if(name == "extmrca.divtime") "divtime" else name, "_rel_error.pdf"), width = 5,height = 7.25)
    } else {
      p_rel_error = p_rel_error +
        theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), legend.position = "none") +
        scale_colour_manual(values=cbPalette) + scale_shape_manual(values = c(16,17))
      
      if(is.null(plotfolder)) show(p_rel_error)
      else ggsave(paste0(plotfolder, "/", if(name == "extmrca.divtime") "divtime" else name, "_rel_error.pdf"),width = 5,height = 5)
    }
    
    p_coverage <- ggplot(get(paste0(name, "_coverage")), mapping = aes(x=Age, y=Coverage, colour=factor(Char_length), shape=factor(Type))) +
      xlim("*Interval ages", "*Symmetric ages", "Correct age", "Median age", "Random age") +
      ylim(0,1) + xlab(plotnames[name]) + ylab("Coverage") + geom_point(size = 3, position = position_dodge(width=0.75)) +
      theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1), text = element_text(size = 15), legend.position = "none") +
      scale_colour_manual(values=cbPalette) + scale_shape_manual(values = c(16, 17))
    
    if(is.null(plotfolder)) show(p_coverage)
    else ggsave(paste0(plotfolder, "/", if(name == "extmrca.divtime") "divtime" else name, "_coverage.pdf"),width = 5,height = 5)
  }
}
