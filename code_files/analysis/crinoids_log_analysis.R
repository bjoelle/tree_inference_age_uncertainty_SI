FBD_analysis = function(treefile, output_folder, highfs = T, relaxed_clock = F, save_folder = "../results/") {
  if(!file.exists(paste0(save_folder,basename(output_folder),"_summary.RData"))) {
    tmpfile = paste0(basename(output_folder),".RData")
    repeat {
      partial1 = FBD_allnodes_divtimes(treefile, output_folder, tmp_file = tmpfile)
      if(class(partial1) != "try-error") break
    } 
    partial2 = FBD_log_analysis(treefile, output_folder, highfs, relaxed_clock)
    partial = c(partial1, partial2)
    psummary = FBD_summary(partial)
    
    names = c("divtime", "diversificationRateFBD","turnoverFBD", "samplingProportionFBD", "SACountFBD", "clockRate")
    full = fsummary = list()
    for (nm in names) {
      full[[nm]] = partial[which(grepl(nm, names(partial)))]
      fsummary[[nm]] = psummary[which(grepl(nm, names(psummary)))]
    }
    full$ESS = partial$ESS
    fsummary$ESS = psummary$ESS
    
    assign(paste0(basename(output_folder),"_analysis"), full)
    save(list = c(paste0(basename(output_folder),"_analysis")), file = paste0(save_folder,basename(output_folder),"_analysis.RData"))
    assign(paste0(basename(output_folder),"_summary"), fsummary)
    save(list = c(paste0(basename(output_folder),"_summary")), file = paste0(save_folder,basename(output_folder),"_summary.RData"))
    file.remove(tmpfile)
  }
}

FBD_summary = function(analysis, ignore_NA = F) {
  res = list()
  if(!is.null(analysis$ESS)) {
    res$ESS = list(npb = analysis$ESS$npb, lpb = analysis$ESS$lpb)
    analysis$ESS = NULL
  }
  if(!is.null(analysis$CLOCK)) {
    res$CLOCK = list(npb = analysis$CLOCK$npb, lpb = analysis$CLOCK$lpb)
    analysis$CLOCK = NULL
  }
  
  for(i in 1:length(analysis)) {
    if(length(analysis[[i]]) > 1) res[[names(analysis)[i]]] = mean(analysis[[i]], na.rm = ignore_NA)
    else res[[names(analysis)[i]]] = analysis[[i]]
  }
  res
}

FBD_log_analysis = function(treefile, output_folder, highfs = F) {
  library(ape)
  library(coda)
  library(Metrics)
  
  spec_rate = 0.06
  ext_rate = 0.045
  sampl_rate = if(highfs) 0.1 else 0.03
  clock_rate = 0.033
  
  load(treefile)
  
  names = c("correct_ages", "median_ages", "interval_ages", "random_ages", "norm_intvl_ages")
  FBDp_names = c("diversificationRateFBD","turnoverFBD", "samplingProportionFBD", "SACountFBD", "clockRate")
  
  trueFBDp = c(spec_rate - ext_rate, ext_rate / spec_rate, sampl_rate / (ext_rate + sampl_rate), 0, clock_rate)
  names(trueFBDp) = FBDp_names
  
  analysis = list()
  ESS = list(npb = 0, lpb = c(), whys = list())
  CLOCK = list(npb = 0, lbp = c(), whens = list())
  
  for(i in 1:length(trees)) {
    print(i)
    logs = list()
    for(x in 1:length(names)) {
      #print(names[[x]])
      logs[[x]] = read.table(file.path(output_folder, paste0("FBD_", names[[x]],"_",i,".log")),header = T,comment.char = "#")
      n = nrow(logs[[x]])
      logs[[x]] = logs[[x]][(round(n/10)+1):n,]
    }
    
    why = c()
    when = c()
    tree = samp_trees[[i]]
    
    #true substitution rate & SA count
    trueFBDp['SACountFBD'] = 0
    for(tip in 1:length(tree$tip.label)) {
      if(tree$edge.length[which(tree$edge[,2] == tip)] == 0) trueFBDp['SACountFBD'] = trueFBDp['SACountFBD'] +1
    }
    
    for(x in 1:length(names)) {
      ess = coda::effectiveSize(as.mcmc(logs[[x]]$posterior))
      #if(ess < 200) browser()
      for(y in FBDp_names) {
        #error on FBD parameters
        if(y == 'SACountFBD') {
          analysis[[paste0(names[x],".",y,".relative_error")]] = c(analysis[[paste0(names[x],".",y,".relative_error")]], abs(median(logs[[x]][,y]) - trueFBDp[y])/length(tree$tip.label))
        } else {
          analysis[[paste0(names[x],".",y,".relative_error")]] = c(analysis[[paste0(names[x],".",y,".relative_error")]], abs(median(logs[[x]][,y]) - trueFBDp[y])/trueFBDp[y])            
          analysis[[paste0(names[x],".",y,".median")]] = c(analysis[[paste0(names[x],".",y,".median")]], median(logs[[x]][,y]))
        }
        HPD = coda::HPDinterval(as.mcmc(logs[[x]][,y]))
        analysis[[paste0(names[x],".",y,".coverage")]] = c(analysis[[paste0(names[x],".",y,".coverage")]], (trueFBDp[y] >= HPD[1]) && (trueFBDp[y] <= HPD[2]))
        analysis[[paste0(names[x],".",y,".HPD_width")]] = c(analysis[[paste0(names[x],".",y,".HPD_width")]], (HPD[2] - HPD[1])/trueFBDp[y] )
        
        if(y == 'clockRate') {
          if(is.infinite(median(logs[[x]][,y]))) when = c(when, paste0(names[x],".",y))
        } else {
          #ESS checks on FBD pars
          ess = coda::effectiveSize(as.mcmc(logs[[x]][,y]))
          if(is.na(ess) || ess<200) why = c(why,paste0(names[x],".",y))
        }
      }
    }
    
    #ESS check summary
    if(length(why) > 0) {
      ESS$npb = ESS$npb + 1
      ESS$lpb = c(ESS$lpb,i)
      ESS$whys[[ESS$npb]] = why
    }
    if(length(when) > 0) {
      CLOCK$npb = CLOCK$npb + 1
      CLOCK$lpb = c(CLOCK$lpb, i)
      CLOCK$whens[[CLOCK$npb]] = when
    }
    
  }
  
  for(x in 1:length(names)) {
    for(y in c("diversificationRateFBD","turnoverFBD", "samplingProportionFBD")) {
      analysis[[paste0(names[x],".",y,".rmse")]] = Metrics::rmse(analysis[[paste0(names[x],".",y,".median")]], rep(trueFBDp[y], length(trees)))
    }
  }
  
  analysis = lapply(analysis, function(x) {names(x) = NULL; x})
  analysis$ESS = ESS
  analysis$CLOCK = CLOCK
  analysis
}

