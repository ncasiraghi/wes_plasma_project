library( data.table )
library( stringr )

filter <- data.frame(max.af.ctrl = rep(0.01,15),
                     min.af.case = c(rep(0.01,3),rep(0.05,12)),
                     min.cov = c(rep(10,6),rep(20,3),rep(30,3),rep(50,3)),
                     min.tc = c(0,0.1,0.2),
                     stringsAsFactors = F )

setwd("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Strelka2_mutation_load_tables/")

toolName <- "Strelka2"

for(indx in seq( 1 ,nrow(filter))){
  
  min.af.case = filter$min.af.case[indx]
  min.cov.case = filter$min.cov[indx] 
  
  max.af.ctrl = filter$max.af.ctrl[indx]
  min.cov.ctrl = filter$min.cov[indx]
  
  min.tc <- filter$min.tc[indx]
  
  label_out <- paste0("afn_",max.af.ctrl,"_aft_",min.af.case,"_mincov_",min.cov.case,"_mintc_",min.tc)
  cat(label_out,"\n")
  
  # Tissue HaloPlex Strelka2 
  halo.calls = readLines(con = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Strelka2/Tissue_HaloPlex/calls_strelka2.txt")
  
  halo.calls = data.frame(id.patient = str_split_fixed(halo.calls, "/", 13)[,12], calls = halo.calls, stringsAsFactors = F)
  
  # Tissue SureSelect Strelka2 
  sure.calls = readLines(con = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Strelka2/Tissue_SureSelect/calls_strelka2.txt")
  
  sure.calls = data.frame(id.patient = str_split_fixed(sure.calls, "/", 13)[,12], calls = sure.calls, stringsAsFactors = F)
  
  # Plasma Strelka2
  calls = readLines(con = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Strelka2/Plasma/calls_strelka2.txt")
  
  plasma.calls = data.frame(id.case = str_split_fixed(calls, "/", 13)[,12], calls = calls, stringsAsFactors = F)
  plasma.calls$id.patient = gsub(plasma.calls$id.case, pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = "")
  
  # Annotations
  anns = read.delim(file = "/elaborazioni/sharedCO/PCF_Project/Analysis_Final/Tables_FREEZE/Annotations.txt",as.is = T,stringsAsFactors = F)
  anns = anns[which(anns$class %in% c("CRPC-Adeno","CRPC-NE","HNPC")),]
  anns = anns[-which(anns$sample %in% c("PM1388","PM1388_Z3_1_Case")),]
  anns$sample = gsub(anns$sample, pattern = "_HALO|-HALO", replacement = "") # clean name for PM12 
  anns.plasma = unique(anns$sample[which(anns$type=="Plasma")])
  anns.tissue = unique(anns$sample[which(anns$type=="Tissue")])
  
  # check and filter by Annotations file 
  plasma.calls = plasma.calls[which(plasma.calls$id.patient %in% unique(gsub(anns.plasma,pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = ""))),]
  
  halo.calls$id.case = gsub(halo.calls$id.patient, pattern = "_HALO|-HALO", replacement = "")
  halo.calls = halo.calls[which(halo.calls$id.case %in% anns.tissue),]
  halo.calls$id.patient = str_split_fixed(halo.calls$id.patient, "_", n = 2)[,1]
  
  sure.calls$id.case = sure.calls$id.patient
  sure.calls = sure.calls[which(sure.calls$id.case %in% anns.tissue),]
  sure.calls$id.patient = str_split_fixed(sure.calls$id.patient, "_", n = 2)[,1]
  
  # < check plasma samples and calls>
  length(anns.plasma)
  nrow(plasma.calls)
  length(anns.plasma) == nrow(plasma.calls)
  setequal(unique(gsub(anns.plasma,pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = "")),plasma.calls$id.patient)
  
  # < check tissue samples and calls>
  length(anns.tissue)
  nrow(halo.calls)+nrow(sure.calls)
  length(anns.tissue) == nrow(halo.calls)+nrow(sure.calls)
  setequal(anns.tissue,c(halo.calls$id.case,sure.calls$id.case))
  
  # global admixture table
  gadm = read.delim("/elaborazioni/sharedCO/PCF_Project/Analysis_Final/Tables_FREEZE/globalAdmTable.txt",stringsAsFactors = F,check.names = F,as.is = T)
  gadm$sample = gsub(gadm$sample, pattern = "_HALO|-HALO", replacement = "") # clean name for PM12 
  gadm <- gadm[which(gadm$sample %in% c(anns.plasma, anns.tissue)), ]
  gadm$purity = 1-gadm$adm
  gadm$purity[which(is.na(gadm$purity))] = 0
  
  plasma.calls <- plasma.calls[which(plasma.calls$id.case %in% gsub( gadm$sample[which( gadm$purity >= min.tc )], pattern = "-1st",replacement = "")),]
  halo.calls <- halo.calls[which(halo.calls$id.case %in% gadm$sample[which( gadm$purity >= min.tc )]),]
  sure.calls <- sure.calls[which(sure.calls$id.case %in% gadm$sample[which( gadm$purity >= min.tc )]),]
  
  # get number of all and nonsyn-only snvs for each sample
  
  get_snvs_count <- function(df, min.af.case, min.cov.case, max.af.ctrl, min.cov.ctrl){
    df$n.snvs <- NA
    df$n.snvs_nonsyn <- NA
    for(i in 1:nrow(df)){
      x <- fread(input = df$calls[i],header = T,stringsAsFactors = F,data.table = F)
      x <- x[which(x$af_case > min.af.case & x$cov_case > min.cov.case),]
      x <- x[which(x$af_control < max.af.ctrl & x$cov_control > min.cov.ctrl),]
      x.nonsyn <- x[which(x$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation")),]
      df$n.snvs[i] <- length(unique(x$group))
      df$n.snvs_nonsyn[i] <- length(unique(x.nonsyn$group))
    }
    return( df )
  }
  
  plasma.calls <- get_snvs_count(df = plasma.calls, min.af.case=min.af.case,min.cov.case=min.cov.case,max.af.ctrl=max.af.ctrl,min.cov.ctrl=min.cov.ctrl)
  halo.calls <- get_snvs_count(df = halo.calls, min.af.case=min.af.case,min.cov.case=min.cov.case,max.af.ctrl=max.af.ctrl,min.cov.ctrl=min.cov.ctrl)
  sure.calls <- get_snvs_count(df = sure.calls, min.af.case=min.af.case,min.cov.case=min.cov.case,max.af.ctrl=max.af.ctrl,min.cov.ctrl=min.cov.ctrl)
  
  # add kit and type
  plasma.calls$Kit <- "NimbleGen"
  plasma.calls$Type <- "Plasma"
  halo.calls$Kit <- "Haloplex"
  halo.calls$Type <- "Tissue"
  sure.calls$Kit <- "SureSelect"
  sure.calls$Type <- "Tissue"
  
  # add annotations
  anns$id.case = NA
  anns$id.case = gsub(anns$sample,pattern = "-1st",replacement = "")
  
  final = rbind(plasma.calls,
                halo.calls,
                sure.calls)
  
  final = merge(x = final,y = anns,by = "id.case",all.x = T)
  
  # < save outs >
  mutload_name_rdata <- paste0("MutationLoadTable_",toolName,"_",label_out,".RData")
  mutload_name_tsv <- paste0("MutationLoadTable_",toolName,"_",label_out,".tsv")
  write.table(x = final[, c( "id.patient", "sample", "n.snvs", "n.snvs_nonsyn", "Type", "class", "Kit")],file = mutload_name_tsv,quote = F,col.names = T,row.names = F,sep = "\t")
  save(final,file = mutload_name_rdata, compress = T)
  
  # table snvs for plasma and tissue
  
  tabout <- c()
  for(i in 1:nrow(final)){
    x <- fread(input = final$calls[i],header = T,stringsAsFactors = F,data.table = F)
    x <- x[which(x$af_case > min.af.case & x$cov_case > min.cov.case),,drop=F]
    x <- x[which(x$af_control < max.af.ctrl & x$cov_control > min.cov.ctrl),,drop=F]
    if(nrow(x)>0){
      x$id.patient <- final$id.patient[i]
      x$sample <- final$sample[i]
      x$Type <- final$Type[i]
      tabout <- rbind(tabout, x)
    }
  }
  
  # < save outs >
  tabsnvs_name_rdata <- paste0("Table_SNVs_tissue_and_plasma_",toolName,"_",label_out,".RData")
  tabsnvs_name_tsv <- paste0("Table_SNVs_tissue_and_plasma_",toolName,"_",label_out,".tsv")
  
  write.table(x = tabout,file = tabsnvs_name_tsv,quote = F,col.names = T,row.names = F,sep = "\t")
  save(tabout,file = tabsnvs_name_rdata, compress = T)
  
}


