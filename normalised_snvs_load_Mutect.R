wd <- "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_mutation_load_tables"
setwd(wd)

ListMutationLoadTable <- grep(list.files(wd,pattern = "MutationLoadTable_Mutect",full.names = TRUE),pattern = '\\.RData$',value = TRUE)

for( MutationLoadTable in ListMutationLoadTable ){
  
  message( MutationLoadTable )
  
  # load the MutationLoadTable
  load(MutationLoadTable)
  final$ratio_median_all <- NA
  
  final$ratio_median_all[which(final$Type=="Tissue")] <- final$n.snvs[which(final$Type=="Tissue")]/median(final$n.snvs[which(final$Type=="Tissue")])
  
  final$ratio_median_all[which(final$Type=="Plasma")] <- final$n.snvs[which(final$Type=="Plasma")]/median(final$n.snvs[which(final$Type=="Plasma")])
  
  save(final,file = gsub(basename(MutationLoadTable),pattern = "MutationLoadTable_Mutect",replacement = "ratio_mutload_medianload_Mutect"),compress = T)
  
}
