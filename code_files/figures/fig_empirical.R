brachiopods_boxplot = function(outfolder, plotfile = NULL) {
  library(ggplot2)
  estnames = c("clockRate", "diversificationRateFBD", "turnoverFBD", "samplingProportionFBD")
  pnames = c("Clock rate", "Diversification rate", "Turnover", "Sampling proportion", "Root age")
  names = c("median_ages", "random_ages", "interval_ages")
  
  logs = list()
  for(x in names) {
    logs[[x]] = read.table(file.path(outfolder, paste0("brachiopods_FBD_", x,".log")),header = T,comment.char = "#")
    n = nrow(logs[[x]])
    logs[[x]] = logs[[x]][(round(n/10)+1):n,]
  }
  
  cbPalette <- c("#56B4E9", "#009E73", "#E69F00", "#CC79A7")
  p = list()
  for(i in 1:length(estnames)) {
    e = estnames[i]
    df = data.frame(values = c(), ds = c())
    for(x in names) {
      df = rbind(df, data.frame(values = logs[[x]][[e]], ds = x))
    }
    px = ggplot(df, aes(x = ds, y = values, fill = ds)) + geom_boxplot(notch = TRUE) + theme_bw() +
      scale_x_discrete(labels = c("Median ages", "Random ages", "Interval ages")) + xlab("") + ylab("") + 
      ggtitle(pnames[i]) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values=cbPalette)
    p[[i]] = px
  }
  
  offsets = c(365.55,367.762481269822)
  df = data.frame(values = c(), ds = c())
  for(i in 1:2) {
    x = names[i]
    df = rbind(df, data.frame(values = logs[[x]][["TreeHeight"]] + offsets[i], ds = x))
  }
  x = "interval_ages"
  df = rbind(df, data.frame(values = logs[[x]][["TreeHeight"]] + logs[[x]][["offset"]] , ds = x))
  p[[length(pnames)]] = ggplot(df, aes(x = ds, y = values, fill = ds)) + geom_boxplot(notch = TRUE) + theme_bw() +
    scale_x_discrete(labels = c("Median ages", "Random ages", "Interval ages")) + xlab("") + ylab("") + 
    ggtitle(pnames[length(pnames)]) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values=cbPalette)
  
  p = gridExtra::grid.arrange(grobs = p, ncol = 3, nrow = 2)
  if(!is.null(plotfile)) ggsave(plotfile, plot = p,width = 10, height = 10)
  else print(p)
}