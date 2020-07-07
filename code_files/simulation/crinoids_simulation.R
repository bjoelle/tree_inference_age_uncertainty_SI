spec_rate = 0.06
ext_rate = 0.045

highfs = TRUE
clock_prior = TRUE #Lognormal Clock Prior with incorrect mean 
correct_clock_prior = FALSE #Lognormal Clock Prior with correct mean (0.033)

if(clock_prior && correct_clock_prior) clock_prior = FALSE

origin_time = 486 #start of ordovician (486 Ma)
simulated_age = 130 #to end of devonian (356 Ma)
ntrees = 50

sampl_rate = if(highfs) 0.1 else 0.03
min_fossils = if(highfs) 80 else 20
max_fossils = if(highfs) 120 else 40
min_span = if(highfs) 50 else 30 

character.length = 300 #30, 300 or 3000
prop_binary = 1.0 
nchars = if(prop_binary == 1.0) 2 else 3 
clock_rate = 0.033 

xml_template = paste0("templates/crinoids_FBD_template_", 
                      if(clock_prior) "clockprior" else if(correct_clock_prior) "correctclockprior" else "nopriors", 
                      ".xml")
name = paste0(if(highfs) "highfs_" else "lowfs_", 
              if(clock_prior) "clockprior_" else if(correct_clock_prior) "correctclockprior_" else "nopriors_", 
              "charlength", character.length)

out_dir = paste0("../", if(highfs) "highfs_" else "lowfs_", "charlength", character.length, "/", name)
save.file = paste0("../crinoid_", if(highfs) "highfs_" else "lowfs_", "charlength", character.length, "/", name, "_seed459.RData")

pbdbf = "../datasets/pbdb_crinoids.RData"
mean_pbdb_span = 12

#source("auxf.R")
if(!dir.exists(out_dir)) dir.create(out_dir)
set.seed(459)

# simulating trees and fossils with parameters
trees = list()
fossils = list()
samp_trees = list()
while(length(trees) < ntrees) {
  nsim = ntrees - length(trees)  
  trees = c(trees, TreeSim::sim.bd.age(simulated_age, nsim, spec_rate, ext_rate,complete = T))
  for (i in ntrees:(ntrees-nsim+1)) {
    if(class(trees[[i]]) != "phylo") { 
      #if trees[[i]] is numeric and not a tree delete from trees/fossils
      trees = trees[-i] 
      if(i <= length(fossils)) fossils = fossils[-i]
      next
    }
    fossils[[i]] = FossilSim::sim.fossils.poisson(sampl_rate, tree = trees[[i]])
    if(length(fossils[[i]]$edge) < min_fossils || length(fossils[[i]]$edge) > max_fossils) { 
      #if fossils[[i]] does not meet the min/max requirments for number of fossils delete from trees/fossils
      fossils = fossils[-i]
      trees = trees[-i]
    }
    else {
      span = max(fossils[[i]]$hmax) - min(fossils[[i]]$hmin)
      if(span < min_span) { 
        #if fossils[[i]] has sampling ages spanning less than min_span delete from trees/fossils
        fossils = fossils[-i]
        trees = trees[-i]
      }
    }
  }
}

library(xml2)
library(magrittr)
sequences = list()
offsets = list()
for (i in 1:ntrees) {
  # adding uncertainty to fossil ages
  fossils[[i]]$h = (fossils[[i]]$hmin + fossils[[i]]$hmax)/2 + origin_time - simulated_age 
  #average min/max sampling time and shift into appropraite timescale
  fossils[[i]] = fossils[[i]][order(fossils[[i]]$edge, -fossils[[i]]$h),] #sort fossils by edge and descending h 
  intervals = get_intervals_from_PBDB(fossils[[i]])
  max.ages = intervals$max
  min.ages = intervals$min
  
  #ftree = combined.tree.with.fossils(trees[[i]],fossils[[i]])
  ftree = FossilSim::SAtree.from.fossils(trees[[i]],fossils[[i]]) 
  #created combined SA (sampled ancestor) tree created from trees[[i]] and fossils[[i]]
  tree = FossilSim::sampled.tree.from.combined(ftree, rho = 0) #removes unsampled lineages from a combined tree
  samp_trees[[i]] = tree
  ages = ape::node.depth.edgelength(tree)
  ages = max(ages) - ages
  offsets[[i]] = list(correct_ages = min(fossils[[i]]$h)) 
  
  # simulating sequences on the trees
  sequences[[i]] = sim_character_sequences_varying_rate(samp_trees[[i]], clock_rate, character.length)
  
  # create xml for certain and uncertain runs
  template = read_xml(xml_template)
  dta = xml_find_first(template,"data")
  tax_ref = rep(F,length(tree$tip.label))
  text = text_med = text_rand = ""
  amed = arand = c()
  
  xml_add_child(dta, as_xml_document(list(userDataType = structure(list(), id="StandardData", spec="beast.evolution.datatype.StandardData", 
                                                                   nrOfStates=as.character(nchars)))))
  # handling sequences & tip dates
  for(j in 1:length(tree$tip.label)) {
    tip_id = as.character(tree$tip.label[j])
    
    xml_add_child(dta, as_xml_document(list(sequence = structure(list(), id = paste0("seq_", tip_id), taxon = tip_id, 
                                                                 totalcount=as.character(nchars), 
                                                                 value = sequences[[i]][tip_id]))))
    if(j > 1) text = paste0(text, ",")
    
    text = paste0(text, "\n", tip_id, "=", round(ages[j],9))
    amed = c(amed, (max.ages[j]+min.ages[j])/2) #median age between min/max
    arand = c(arand, runif(1,min.ages[j],max.ages[j])) #random age in min/max range
  }
  
  offset_med = min(amed)
  offset_rand = min(arand)
  for(j in 1:length(tree$tip.label)) {
    if(j > 1) {
      text_med = paste0(text_med, ",")
      text_rand = paste0(text_rand, ",")
    }
    tip_id = as.character(tree$tip.label[j])
    text_rand = paste0(text_rand, "\n", tip_id, "=", arand[j] - offset_rand)
    text_med = paste0(text_med, "\n", tip_id, "=", amed[j] - offset_med)
  }
  offsets[[i]]$median_ages = offset_med
  offsets[[i]]$random_ages = offset_rand
  
  xml_replace(xml_find_first(template,"data"),dta)
  
  # different tip dates & names for median & correct ages
  xml_find_first(template,"run/state/tree/trait") %>% xml_set_text(text_med)
  xml_find_first(template,"run/distribution/distribution/distribution/treeWOffset") %>% xml_set_attr("offset", offsets[[i]]$median_ages) 
  xml_find_first(template,"run/logger") %>% xml_set_attr("fileName", paste0("FBD_median_ages_",i,".log"))
  xml_find_first(template, "run/logger[@mode='tree']") %>% xml_set_attr("fileName", paste0("FBD_median_ages_",i,".trees"))
  write_xml(template,file.path(out_dir,paste0("FBD_median_ages_",i,".xml")))
  
  xml_find_first(template,"run/state/tree/trait") %>% xml_set_text(text_rand)
  xml_find_first(template,"run/distribution/distribution/distribution/treeWOffset") %>% xml_set_attr("offset", offsets[[i]]$random_ages) 
  xml_find_first(template,"run/logger") %>% xml_set_attr("fileName", paste0("FBD_random_ages_",i,".log"))
  xml_find_first(template, "run/logger[@mode='tree']") %>% xml_set_attr("fileName", paste0("FBD_random_ages_",i,".trees"))
  write_xml(template,file.path(out_dir,paste0("FBD_random_ages_",i,".xml")))
  
  xml_find_first(template,"run/state/tree/trait") %>% xml_set_text(text)
  xml_find_first(template,"run/distribution/distribution/distribution/treeWOffset") %>% xml_set_attr("offset", offsets[[i]]$correct_ages) 
  xml_find_first(template,"run/logger") %>% xml_set_attr("fileName", paste0("FBD_correct_ages_",i,".log"))
  xml_find_first(template, "run/logger[@mode='tree']") %>% xml_set_attr("fileName", paste0("FBD_correct_ages_",i,".trees"))
  write_xml(template,file.path(out_dir,paste0("FBD_correct_ages_",i,".xml")))
  
  # adding uncertain ages
  operator = as_xml_document(list(operator = structure(list(), spec="SampledNodeDateRandomWalker", 
                                                       windowSize="10.", tree="@Tree.t:26", treeWOffset="@treeWoffset", weight="5.")))
  operator_n = as_xml_document(list(operator = structure(list(), spec="SampledNodeDateRandomWalker", 
                                                         windowSize="10.", tree="@Tree.t:26", treeWOffset="@treeWoffset", weight="5.")))
  taxset = as_xml_document(list(taxonset = structure(list(), id = "fossilTaxa", spec="TaxonSet")))
  for(j in 1:length(tree$tip.label)) {
    
    tip_id = as.character(tree$tip.label[j])
    if(!tax_ref[j]) xml_add_child(taxset, as_xml_document(list(taxon = structure(list(), id=tip_id, spec="Taxon"))))
    else xml_add_child(taxset, as_xml_document(list(taxon = structure(list(), idref=tip_id))))
    
    xml_add_child(operator, as_xml_document(list(samplingDates = structure(list(),
                                                                           id = paste0("samplingDate.",tip_id), spec="beast.evolution.tree.SamplingDate", taxon = paste0("@",tip_id),
                                                                           lower=as.character(min.ages[j]), upper=as.character(max.ages[j])))))
    xml_add_child(operator_n, as_xml_document(list(samplingDates = structure(list(),
                                                                             id = paste0("samplingDate.",tip_id), spec="beast.evolution.tree.SamplingDate", taxon = paste0("@",tip_id),
                                                                             lower=as.character(max(0, ages[j] + offsets[[i]]$correct_ages - mean_pbdb_span/2)), upper=as.character(ages[j] + offsets[[i]]$correct_ages + min(ages[j] + offsets[[i]]$correct_ages,mean_pbdb_span/2))))))
  }
  
  xml_add_child(operator, taxset)
  xml_find_first(template,"run") %>% xml_add_child(operator)
  
  xml_find_first(template,"run/logger") %>% xml_set_attr("fileName", paste0("FBD_interval_ages_",i,".log"))
  xml_find_first(template, "run/logger[@mode='tree']") %>% xml_set_attr("fileName", paste0("FBD_interval_ages_",i,".trees"))
  write_xml(template,file.path(out_dir,paste0("FBD_interval_ages_",i,".xml")))
  
  xml_add_child(operator_n, taxset)
  xml_replace(xml_find_first(template,"run/operator[@spec='SampledNodeDateRandomWalker']"),operator_n)
  
  xml_find_first(template,"run/logger") %>% xml_set_attr("fileName", paste0("FBD_norm_intvl_ages_",i,".log"))
  xml_find_first(template, "run/logger[@mode='tree']") %>% xml_set_attr("fileName", paste0("FBD_norm_intvl_ages_",i,".trees"))
  write_xml(template,file.path(out_dir,paste0("FBD_norm_intvl_ages_",i,".xml")))
  
  # write nexus files for unconstrained analysis
  ape::write.nexus.data(strsplit(sequences[[i]], split = ""), file = paste0(out_dir, "/unconstrained/rep_", i, ".nex"), format = "standard", gap = "-", missing = "?")  
}

save(trees, fossils, samp_trees, offsets, sequences, file = save.file)
