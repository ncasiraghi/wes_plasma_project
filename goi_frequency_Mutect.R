wd <- "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_mutation_load_tables"
setwd(wd)

MutationLoadTable <- "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_mutation_load_tables/MutationLoadTable_Mutect_afn_0.01_aft_0.05_mincov_30_mintc_0.RData"

# Genes of interest
goi = read.delim(file = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/abemus_plasma_tissues_distances/Genes_List_Plasma_Study_180321.txt",as.is = T,stringsAsFactors = F)
to.exclude = c("FANCA","ERG","MET","BRAF","CYLD","AR")
goi = rbind(goi,c("PMS2",rep(NA,5)))
goi = goi[which(!goi$symbol%in%to.exclude),]

# Annotations
anns = read.delim(file = "/elaborazioni/sharedCO/PCF_Project/Analysis_Final/Tables_FREEZE/Annotations.txt",as.is = T,stringsAsFactors = F)
anns = anns[which(anns$class %in% c("CRPC-Adeno","CRPC-NE","HNPC")),]
anns = anns[-which(anns$sample %in% c("PM1388","PM1388_Z3_1_Case")),]
anns$sample = gsub(anns$sample, pattern = "_HALO|-HALO", replacement = "") # clean name for PM12 
anns.plasma = unique(anns$sample[which(anns$type=="Plasma")])
anns.tissue = unique(anns$sample[which(anns$type=="Tissue")])

hnpc = unique(gsub(anns$sample[which(anns$type=="Plasma" & anns$class=="HNPC")],pattern = "-1st",replacement = ""))
nepc = unique(gsub(anns$sample[which(anns$type=="Plasma" & anns$class=="CRPC-NE")],pattern = "-1st",replacement = ""))
adeno= unique(gsub(anns$sample[which(anns$type=="Plasma" & anns$class=="CRPC-Adeno")],pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = ""))
#adeno = setdiff(adeno,c("PM14","PM185"))

# load the MutationLoadTable
load(MutationLoadTable)
final <- final[which(!final$id.case %in% c("PM14-TP2","PM14-TP3","PM14-TP4","PM185-TP2","PM185-TP3")),]

final$class_plasma_based = NA
final$class_plasma_based[which(final$id.patient %in% hnpc)] <- "HNPC"
final$class_plasma_based[which(final$id.patient %in% nepc)] <- "CRPC-NE"
final$class_plasma_based[which(final$id.patient %in% adeno)]<- "CRPC-Adeno"

# Frequencies of aberration for Gene List
get.final.freq <- function(data_out){
  
  n = 7:ncol(data_out)
  count.miss=apply(X = data_out[,n],MARGIN = 1,function(x) sum(x==1,na.rm=T))
  count.nons=apply(X = data_out[,n],MARGIN = 1,function(x) sum(x==2,na.rm=T))
  count.both=apply(X = data_out[,n],MARGIN = 1,function(x) sum(x==3,na.rm=T))
  count=count.miss+count.nons+count.both
  
  freq = count/length(n)
  freq.miss = count.miss/length(n)
  freq.nons = count.nons/length(n)
  freq.both = count.both/length(n)
  
  data_out$aberr.count = count
  data_out$aberr.freq = freq
  data_out$aberr.count.miss = count.miss
  data_out$aberr.freq.miss = freq.miss
  data_out$aberr.count.nons = count.nons
  data_out$aberr.freq.nons = freq.nons
  data_out$aberr.count.both = count.both
  data_out$aberr.freq.both = freq.both
  
  return(data_out)
  
}

getAberrantGenes = function(dataset, goi){
  out <- goi
  out.adeno <- goi
  out.nepc <- goi
  for(i in 1:nrow(dataset)){
    
    m <- read.delim(dataset$calls[i],as.is = T,stringsAsFactors = F)
    m <- m[which(m$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation")),]
    
    w = as.data.frame(table(m$Hugo_Symbol))
    w.miss = as.data.frame(table(m$Hugo_Symbol[which(m$Variant_Classification == "Missense_Mutation")]))
    w.nons = as.data.frame(table(m$Hugo_Symbol[which(m$Variant_Classification == "Nonsense_Mutation")]))
    
    if(nrow(w.nons)==0 & nrow(w.miss)==0){
      w <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(w) = c("symbol",dataset$id.case[i])
    }
    if(nrow(w.nons)==0 & nrow(w.miss)>0){
      colnames(w) = c("symbol",dataset$id.case[i])
      w[,2] <- 1
    }
    if(nrow(w.miss)==0 & nrow(w.nons)>0){
      colnames(w) = c("symbol",dataset$id.case[i])
      w[,2] <- 2
    }
    if(nrow(w.miss)>0 & nrow(w.nons)>0){
      colnames(w) = c("symbol",dataset$id.case[i])
      colnames(w.miss) = c("symbol",dataset$id.case[i])
      colnames(w.nons) = c("symbol",dataset$id.case[i])
      w.miss[,2] = 1
      w.nons[,2] = 2
      x = merge(x = w.miss,y = w.nons,by = "symbol",all = T,suffixes = c("_miss","_nons"))
      w[,2] = apply(X = x[,2:3],MARGIN = 1,FUN = sum,na.rm=T)
    }
    
    out = merge(x = out,y = w,by = "symbol",all.x = T)
    
    sample.class <- dataset$class_plasma_based[i]
    
    if(sample.class == "CRPC-NE"){
      out.nepc = merge(x = out.nepc,y = w,by = "symbol",all.x = T)
    }
    if(sample.class == "CRPC-Adeno"){
      out.adeno = merge(x = out.adeno,y = w,by = "symbol",all.x = T)
    }
    
  }
  
  out <- get.final.freq(data_out = out)
  out.adeno <- get.final.freq(data_out = out.adeno)
  out.nepc <- get.final.freq(data_out = out.nepc)
  
  return(list(out,out.adeno,out.nepc))
}


# Data to plot for aberration frequency 
plasma <- final[which(final$class_plasma_based %in% c("CRPC-Adeno","CRPC-NE") & final$Type == "Plasma") ,]
res = getAberrantGenes(dataset = plasma, goi = goi)
aberrant.plasma = res[[1]]
row.names(aberrant.plasma) = aberrant.plasma$symbol
aberrant.plasma.adeno = res[[2]]
row.names(aberrant.plasma.adeno) = aberrant.plasma.adeno$symbol
aberrant.plasma.nepc = res[[3]]
row.names(aberrant.plasma.nepc) = aberrant.plasma.nepc$symbol

tissue <- final[which(final$class_plasma_based %in% c("CRPC-Adeno","CRPC-NE") & final$Type == "Tissue") ,]
res = getAberrantGenes(dataset = tissue, goi = goi)
aberrant.tissue = res[[1]]
row.names(aberrant.tissue) = aberrant.tissue$symbol
aberrant.tissue.adeno = res[[2]]
row.names(aberrant.tissue.adeno) = aberrant.tissue.adeno$symbol
aberrant.tissue.nepc = res[[3]]
row.names(aberrant.tissue.nepc) = aberrant.tissue.nepc$symbol

a = aberrant.plasma[,c("symbol","aberr.freq")]
b = aberrant.tissue[,c("symbol","aberr.freq")]
tmp = merge(x = a,y = b,by = "symbol",suffixes = c("_plasma","_tissue"))
tmp=tmp[with(tmp, order(aberr.freq_plasma,aberr.freq_tissue)), ]

tmp = tmp[which(tmp$aberr.freq_plasma>0 | tmp$aberr.freq_tissue>0),]

aberrant.plasma = aberrant.plasma[tmp$symbol, ]

aberrant.plasma.adeno = aberrant.plasma.adeno[row.names(aberrant.plasma),]
aberrant.plasma.nepc = aberrant.plasma.nepc[row.names(aberrant.plasma),]
aberrant.tissue = aberrant.tissue[row.names(aberrant.plasma),]
aberrant.tissue.adeno = aberrant.tissue.adeno[row.names(aberrant.plasma),]
aberrant.tissue.nepc = aberrant.tissue.nepc[row.names(aberrant.plasma),]

save(aberrant.plasma,
     aberrant.plasma.adeno,
     aberrant.plasma.nepc,
     aberrant.tissue,
     aberrant.tissue.adeno,
     aberrant.tissue.nepc,
     file = gsub(basename(MutationLoadTable),pattern = "MutationLoadTable_Mutect",replacement = "goi.aberrant.freq_Mutect"),compress = T)
