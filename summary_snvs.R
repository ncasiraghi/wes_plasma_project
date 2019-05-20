
mutation_load_tables_folder <- "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Strelka2_mutation_load_tables"

setwd(mutation_load_tables_folder)

mlt_list <- grep(list.files(mutation_load_tables_folder,pattern = "MutationLoadTable_",full.names = T),pattern = "\\.RData$",value = T)

for(mlt in mlt_list){
  message(mlt)
  pattern <- gsub(basename(mlt),pattern = "MutationLoadTable_",replacement = "")
  snt <- grep(list.files(mutation_load_tables_folder,pattern = "Table_SNVs_tissue_and_plasma_",full.names = T),pattern = pattern,value = T)
  load(mlt)
  load(snt)
  
  pdf(file = paste0("summary_snvs_",gsub(pattern,pattern = "\\.RData",replacement = ""),".pdf"), h = 219, w = 297, paper='A4r')
  
  layout(matrix(c(1,1,2,2,
                  1,1,2,2,
                  3,3,4,4,
                  3,3,4,4), nrow=4, byrow=TRUE))
  
  par(pty="s")
  
  boxplot(n.snvs ~ class,data = final, varwidth=T, ylab="n.snvs")
  boxplot(n.snvs_nonsyn ~ class, data = final, varwidth=T, ylab="n.snvs_nonsyn")
  boxplot(n.snvs ~ Type, data = final, varwidth=T, ylab="n.snvs")
  boxplot(n.snvs_nonsyn ~ Type, data = final, varwidth=T, ylab="n.snvs_nonsyn")
  
  library( hexbin )
  plasma <- tabout[which(tabout$Type=="Plasma"),]
  tissue <- tabout[which(tabout$Type=="Tissue"),]
  
  bin<-hexbin(x = plasma$af_control,y = plasma$af_case, xbins=50)
  plot(bin,main="all snvs in plasma",xlab="af_control",ylab="af_case") 
  
  bin<-hexbin(x = tissue$af_control,y = tissue$af_case, xbins=50)
  plot(bin,main="all snvs in tissue",xlab="af_control",ylab="af_case") 
  
  dev.off()

}

