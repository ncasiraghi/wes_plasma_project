min.tc = 0
#min.tc = 0.1
#min.tc = 0.2

min.af.case = 0.05

max.af.ctrl = 0.01
min.cov.ctrl = 0

for(min.cov.case in c(10,20,30,50)){
  if(min.tc == 0){
    dir <- file.path("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/", paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case))
    setwd(dir)
  }else{
    dir <- file.path("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/", paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case, "_mintc_", min.tc))
    setwd(dir)
  }
  #setwd("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/")
  
  #dms = list.files(path = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/abemus_plasma_tissues_distances/mintc_0.10_afth_0.05_sf2.4_afthctrl_0.05_dbsnp_pbemAFpar/",pattern = ".RData$",full.names = T)
  dms = list.files(path = dir, pattern = ".RData$",full.names = T)
  
  # annotations
  anns = read.delim(file = "/elaborazioni/sharedCO/PCF_Project/Analysis_Final/Tables_FREEZE/Annotations.txt",as.is = T,stringsAsFactors = F)
  anns.plasma = anns[which(anns$type=="Plasma"),]
  hnpc = unique(gsub(anns.plasma$sample[which(anns.plasma$class == "HNPC")],pattern = "-1st",replacement = ""))
  nepc = unique(gsub(anns.plasma$sample[which(anns.plasma$class == "CRPC-NE")],pattern = "-1st",replacement = ""))
  adeno= unique(gsub(anns.plasma$sample[which(anns.plasma$class == "CRPC-Adeno")],pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = ""))
  adeno = setdiff(adeno,c("PM14","PM185"))
  
  dist.adeno = c()
  dist.nepc = c()
  dist.adeno.plasma_tissue = c()
  dist.nepc.plasma_tissue = c()
  dist.adeno.tissue_plasma = c()
  dist.nepc.tissue_plasma = c()
  
  get.mean.similarity.tissues = function(k){
    out=c()
    if(nrow(k)>0){
      
      for(idx in seq(1,nrow(k),by = 2)){
        data = c(k$similarity[idx],k$similarity[idx+1])
        out = c(out,mean(data))
      }
      
    }
    return(out)
  }
  
  #x=dms[1]
  count.samples = c()
  
  for(x in dms){
    message(basename(x))
    load(x)
    id.patient = gsub(basename(x),pattern = "distances_matrix_|.RData",replacement = "")
    plasma_time_point_ids = c(paste0(id.patient,"-TP1"),paste0(id.patient,"-TP2"),paste0(id.patient,"-TP3"),paste0(id.patient,"-TP4"))
    
    count.plasma = length(which(colnames(m) %in% id.patient))
    count.plasma.timepoint = length(which(colnames(m) %in% plasma_time_point_ids))
    count.tissue = length(colnames(m))-sum(count.plasma,count.plasma.timepoint)
    
    this = data.frame(id.patient=id.patient,
                      count.plasma=count.plasma,
                      count.plasma.timepoint=count.plasma.timepoint,
                      count.tissue=count.tissue,
                      stringsAsFactors = F)
    this$class = NA
    
    head( dd )
    dd$pair = "tissue-tissue"
    dd$pair[ which( dd$a == id.patient) ] <- "plasma-tissue"
    dd$pair[ which( dd$b == id.patient) ] <- "tissue-plasma"
    
    #dd$pair[which(dd$a==id.patient)] = "tissue-plasma"
    #dd$pair[which(dd$b==id.patient)] = "plasma- tissue"
    
    if(id.patient%in%adeno){
      this$class = "adeno"
      dist.adeno = c(dist.adeno,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
      #dist.adeno = c(dist.adeno,dd$similarity[which(dd$pair == "tissue-tissue")])
      dist.adeno.plasma_tissue = c(dist.adeno.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
      dist.adeno.tissue_plasma = c(dist.adeno.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
    }
    if(id.patient%in%nepc){
      this$class = "nepc"
      dist.nepc = c(dist.nepc,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
      #dist.nepc = c(dist.nepc,dd$similarity[which(dd$pair == "tissue-tissue")])
      dist.nepc.plasma_tissue = c(dist.nepc.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
      dist.nepc.tissue_plasma = c(dist.nepc.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
    }
    if(id.patient == "PM14"){
      this$class = "adeno"
      to.remove = unique(c(grep(dd$a,pattern = "-TP2|-TP3|-TP4"),grep(dd$b,pattern = "-TP2|-TP3|-TP4")))
      dd = dd[-to.remove,]
      dd$pair[which(dd$a=="PM14-TP1")] = "tissue-plasma"
      dd$pair[which(dd$b=="PM14-TP1")] = "plasma-tissue"
      dist.adeno = c(dist.adeno,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
      #dist.adeno = c(dist.adeno,dd$similarity[which(dd$pair == "tissue-tissue")])
      dist.adeno.plasma_tissue = c(dist.adeno.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
      dist.adeno.tissue_plasma = c(dist.adeno.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
    }
    if(id.patient == "PM185"){ 
      this$class = "adeno"
      test = which(dd$a == id.patient)
      if(length(test)>0){
        to.remove = unique(c(grep(dd$a,pattern = "-TP1|-TP2|-TP3|-TP4"),grep(dd$b,pattern = "-TP1|-TP2|-TP3|-TP4")))
        dd = dd[-to.remove,]
        dist.adeno = c(dist.adeno,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
        #dist.adeno = c(dist.adeno,dd$similarity[which(dd$pair == "tissue-tissue")])
        dist.adeno.plasma_tissue = c(dist.adeno.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
        dist.adeno.tissue_plasma = c(dist.adeno.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
      } else {
        to.remove = unique(c(grep(dd$a,pattern = "-TP2|-TP3|-TP4"),grep(dd$b,pattern = "-TP2|-TP3|-TP4")))
        dd = dd[-to.remove,]
        dist.adeno = c(dist.adeno,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
        #dist.adeno = c(dist.adeno,dd$similarity[which(dd$pair == "tissue-tissue")])
        dist.adeno.plasma_tissue = c(dist.adeno.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
        dist.adeno.tissue_plasma = c(dist.adeno.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
      }
    }
    count.samples = rbind(count.samples,this)
  }
  
  count.samples = count.samples[which(!is.na(count.samples$class)),]
  count.samples$count.plasma[which(count.samples$id.patient=="PM14")] = 1
  count.samples$count.plasma[which(count.samples$id.patient=="PM185")] = 1
  
  if(min.tc == 0){
    pdf(paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case, "__similarities_distributions_sameaxes.pdf"), width = 8, height = 4)
  }else{
    pdf(paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case, "_mintc_", min.tc, "__similarities_distributions_sameaxes.pdf"), width = 8, height = 4)
  }
  
  
  nepc.color <- c("#682A37")
  nepc.color.rgb = rgb(104,42,55,maxColorValue = 255)
  crpc.color <- c("#BF9C96")
  crpc.color.rgb = rgb(191,156,150,maxColorValue = 255)
  
  par(mfrow=c(1,3),mar=c(4,4,4,1))
  boxplot(dist.adeno,dist.nepc,varwidth = T,names = c("CRPC-Adeno","CRPC-NE"),main="Among Tissues",ylab="SNV Similarity (Non-synonymous)",ylim=c(0,1),outline = F,col=c(crpc.color,nepc.color))
  stripchart(list(dist.adeno,dist.nepc), vertical=T, method="jitter",add=TRUE,pch=19,col=rgb(128,128,128,maxColorValue = 255),jitter = 0.2)
  p = wilcox.test(x = dist.adeno,y = dist.nepc)
  p = formatC(p$p.value,format="e",digits=2)
  text(1.5,1,paste("p = ",p,sep=""))
  
  text(1.5,0.96,paste("n.patients =",paste(length(unique(count.samples$id.patient[which(count.samples$class=="adeno" & count.samples$count.tissue>1)])),
                                           length(unique(count.samples$id.patient[which(count.samples$class=="nepc" & count.samples$count.tissue>1)])),sep="|")))
  
  text(1.5,0.91,paste("n.tissues =",paste(sum(count.samples$count.tissue[which(count.samples$class=="adeno" & count.samples$count.tissue>1)]),
                                          sum(count.samples$count.tissue[which(count.samples$class=="nepc" & count.samples$count.tissue>1)]),sep="|")))
  
  text(1.5,0.86,paste("n.points = ",paste(length(dist.adeno),length(dist.nepc),sep="|")))
  
  
  boxplot(dist.adeno.plasma_tissue,dist.nepc.plasma_tissue,varwidth = T,names = c("CRPC-Adeno","CRPC-NE"),main="Fraction of plasma in tissue",ylab="SNV Similarity (Non-synonymous)",ylim=c(0,1),outline = F) #ylim=c(0,0.3)
  stripchart(list(dist.adeno.plasma_tissue,dist.nepc.plasma_tissue), vertical=T, method="jitter",add=TRUE,pch=19,col=c(crpc.color.rgb,nepc.color.rgb),jitter = 0.2)
  p = wilcox.test(x = dist.adeno.plasma_tissue,dist.nepc.plasma_tissue)
  p = formatC(p$p.value,format="e",digits=2)
  text(1.5,1,paste("p = ",p,sep=""))
  
  text(1.5,1*0.96,paste("n.patients =",paste(length(unique(count.samples$id.patient[which(count.samples$class=="adeno")])),
                                                    length(unique(count.samples$id.patient[which(count.samples$class=="nepc")])),sep="|")))
  
  text(1.5,1*0.91,paste("n.plasma =",paste(sum(count.samples$count.plasma[which(count.samples$class=="adeno")]),
                                                  sum(count.samples$count.plasma[which(count.samples$class=="nepc")]),sep="|")))
  
  text(1.5,1*0.86,paste("n.tissues =",paste(sum(count.samples$count.tissue[which(count.samples$class=="adeno")]),
                                                   sum(count.samples$count.tissue[which(count.samples$class=="nepc")]),sep="|")))
  
  text(1.5,1*0.81,paste("n.points = ",paste(length(dist.adeno.plasma_tissue),length(dist.nepc.plasma_tissue),sep="|")))
  
  boxplot(dist.adeno.tissue_plasma,dist.nepc.tissue_plasma,varwidth = T,names = c("CRPC-Adeno","CRPC-NE"),main="Fraction of tissue in plasma",ylab="SNV Similarity (Non-synonymous)",ylim=c(0,1),outline = F)
  stripchart(list(dist.adeno.tissue_plasma,dist.nepc.tissue_plasma), vertical=T, method="jitter",add=TRUE,pch=19,col=c(crpc.color.rgb,nepc.color.rgb),jitter = 0.2)
  p = wilcox.test(x = dist.adeno.tissue_plasma,dist.nepc.tissue_plasma)
  p = formatC(p$p.value,format="e",digits=2)
  text(1.5,1,paste("p = ",p,sep=""))
  
  text(1.5,0.96,paste("n.patients =",paste(length(unique(count.samples$id.patient[which(count.samples$class=="adeno")])),
                                           length(unique(count.samples$id.patient[which(count.samples$class=="nepc")])),sep="|")))
  
  text(1.5,0.91,paste("n.plasma =",paste(sum(count.samples$count.plasma[which(count.samples$class=="adeno")]),
                                         sum(count.samples$count.plasma[which(count.samples$class=="nepc")]),sep="|")))
  
  text(1.5,0.86,paste("n.tissues =",paste(sum(count.samples$count.tissue[which(count.samples$class=="adeno")]),
                                          sum(count.samples$count.tissue[which(count.samples$class=="nepc")]),sep="|")))
  
  text(1.5,0.81,paste("n.points = ",paste(length(dist.adeno.tissue_plasma),length(dist.nepc.tissue_plasma),sep="|")))
  
  dev.off()
}


min.af.case = 0.01
min.cov.case = 10


if(min.tc == 0){
  dir <- file.path("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/", paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case))
  setwd(dir)
}else{
  dir <- file.path("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/", paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case, "_mintc_", min.tc))
  setwd(dir)
}
#setwd("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/")

#dms = list.files(path = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/abemus_plasma_tissues_distances/mintc_0.10_afth_0.05_sf2.4_afthctrl_0.05_dbsnp_pbemAFpar/",pattern = ".RData$",full.names = T)
dms = list.files(path = dir, pattern = ".RData$",full.names = T)

# annotations
anns = read.delim(file = "/elaborazioni/sharedCO/PCF_Project/Analysis_Final/Tables_FREEZE/Annotations.txt",as.is = T,stringsAsFactors = F)
anns.plasma = anns[which(anns$type=="Plasma"),]
hnpc = unique(gsub(anns.plasma$sample[which(anns.plasma$class == "HNPC")],pattern = "-1st",replacement = ""))
nepc = unique(gsub(anns.plasma$sample[which(anns.plasma$class == "CRPC-NE")],pattern = "-1st",replacement = ""))
adeno= unique(gsub(anns.plasma$sample[which(anns.plasma$class == "CRPC-Adeno")],pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = ""))
adeno = setdiff(adeno,c("PM14","PM185"))

dist.adeno = c()
dist.nepc = c()
dist.adeno.plasma_tissue = c()
dist.nepc.plasma_tissue = c()
dist.adeno.tissue_plasma = c()
dist.nepc.tissue_plasma = c()

get.mean.similarity.tissues = function(k){
  out=c()
  if(nrow(k)>0){
    
    for(idx in seq(1,nrow(k),by = 2)){
      data = c(k$similarity[idx],k$similarity[idx+1])
      out = c(out,mean(data))
    }
    
  }
  return(out)
}

#x=dms[1]
count.samples = c()

for(x in dms){
  message(basename(x))
  load(x)
  id.patient = gsub(basename(x),pattern = "distances_matrix_|.RData",replacement = "")
  plasma_time_point_ids = c(paste0(id.patient,"-TP1"),paste0(id.patient,"-TP2"),paste0(id.patient,"-TP3"),paste0(id.patient,"-TP4"))
  
  count.plasma = length(which(colnames(m) %in% id.patient))
  count.plasma.timepoint = length(which(colnames(m) %in% plasma_time_point_ids))
  count.tissue = length(colnames(m))-sum(count.plasma,count.plasma.timepoint)
  
  this = data.frame(id.patient=id.patient,
                    count.plasma=count.plasma,
                    count.plasma.timepoint=count.plasma.timepoint,
                    count.tissue=count.tissue,
                    stringsAsFactors = F)
  this$class = NA
  
  head( dd )
  dd$pair = "tissue-tissue"
  dd$pair[ which( dd$a == id.patient) ] <- "plasma-tissue"
  dd$pair[ which( dd$b == id.patient) ] <- "tissue-plasma"
  
  #dd$pair[which(dd$a==id.patient)] = "tissue-plasma"
  #dd$pair[which(dd$b==id.patient)] = "plasma-tissue"
  
  if(id.patient%in%adeno){
    this$class = "adeno"
    dist.adeno = c(dist.adeno,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
    #dist.adeno = c(dist.adeno,dd$similarity[which(dd$pair == "tissue-tissue")])
    dist.adeno.plasma_tissue = c(dist.adeno.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
    dist.adeno.tissue_plasma = c(dist.adeno.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
  }
  if(id.patient%in%nepc){
    this$class = "nepc"
    dist.nepc = c(dist.nepc,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
    #dist.nepc = c(dist.nepc,dd$similarity[which(dd$pair == "tissue-tissue")])
    dist.nepc.plasma_tissue = c(dist.nepc.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
    dist.nepc.tissue_plasma = c(dist.nepc.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
  }
  if(id.patient == "PM14"){
    this$class = "adeno"
    to.remove = unique(c(grep(dd$a,pattern = "-TP2|-TP3|-TP4"),grep(dd$b,pattern = "-TP2|-TP3|-TP4")))
    dd = dd[-to.remove,]
    dd$pair[which(dd$a=="PM14-TP1")] = "tissue-plasma"
    dd$pair[which(dd$b=="PM14-TP1")] = "plasma-tissue"
    dist.adeno = c(dist.adeno,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
    #dist.adeno = c(dist.adeno,dd$similarity[which(dd$pair == "tissue-tissue")])
    dist.adeno.plasma_tissue = c(dist.adeno.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
    dist.adeno.tissue_plasma = c(dist.adeno.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
  }
  if(id.patient == "PM185"){
    this$class = "adeno"
    test = which(dd$a == id.patient)
    if(length(test)>0){
      to.remove = unique(c(grep(dd$a,pattern = "-TP1|-TP2|-TP3|-TP4"),grep(dd$b,pattern = "-TP1|-TP2|-TP3|-TP4")))
      dd = dd[-to.remove,]
      dist.adeno = c(dist.adeno,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
      #dist.adeno = c(dist.adeno,dd$similarity[which(dd$pair == "tissue-tissue")])
      dist.adeno.plasma_tissue = c(dist.adeno.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
      dist.adeno.tissue_plasma = c(dist.adeno.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
    } else {
      to.remove = unique(c(grep(dd$a,pattern = "-TP2|-TP3|-TP4"),grep(dd$b,pattern = "-TP2|-TP3|-TP4")))
      dd = dd[-to.remove,]
      dist.adeno = c(dist.adeno,get.mean.similarity.tissues(k = dd[which(dd$pair == "tissue-tissue"),,drop=F]))
      #dist.adeno = c(dist.adeno,dd$similarity[which(dd$pair == "tissue-tissue")])
      dist.adeno.plasma_tissue = c(dist.adeno.plasma_tissue,dd$similarity[which(dd$pair == "plasma-tissue")])
      dist.adeno.tissue_plasma = c(dist.adeno.tissue_plasma,dd$similarity[which(dd$pair == "tissue-plasma")])
    }
  }
  count.samples = rbind(count.samples,this)
}

count.samples = count.samples[which(!is.na(count.samples$class)),]
count.samples$count.plasma[which(count.samples$id.patient=="PM14")] = 1
count.samples$count.plasma[which(count.samples$id.patient=="PM185")] = 1

if(min.tc == 0){
    pdf(paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case, "__similarities_distributions_sameaxes.pdf"), width = 8, height = 4)
  }else{
    pdf(paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case, "_mintc_", min.tc, "__similarities_distributions_sameaxes.pdf"), width = 8, height = 4)
  }

nepc.color <- c("#682A37")
nepc.color.rgb = rgb(104,42,55,maxColorValue = 255)
crpc.color <- c("#BF9C96")
crpc.color.rgb = rgb(191,156,150,maxColorValue = 255)

par(mfrow=c(1,3),mar=c(4,4,4,1))
boxplot(dist.adeno,dist.nepc,varwidth = T,names = c("CRPC-Adeno","CRPC-NE"),main="Among Tissues",ylab="SNV Similarity (Non-synonymous)",ylim=c(0,1),outline = F,col=c(crpc.color,nepc.color))
stripchart(list(dist.adeno,dist.nepc), vertical=T, method="jitter",add=TRUE,pch=19,col=rgb(128,128,128,maxColorValue = 255),jitter = 0.2)
p = wilcox.test(x = dist.adeno,y = dist.nepc)
p = formatC(p$p.value,format="e",digits=2)
text(1.5,1,paste("p = ",p,sep=""))

text(1.5,0.96,paste("n.patients =",paste(length(unique(count.samples$id.patient[which(count.samples$class=="adeno" & count.samples$count.tissue>1)])),
                                         length(unique(count.samples$id.patient[which(count.samples$class=="nepc" & count.samples$count.tissue>1)])),sep="|")))

text(1.5,0.91,paste("n.tissues =",paste(sum(count.samples$count.tissue[which(count.samples$class=="adeno" & count.samples$count.tissue>1)]),
                                        sum(count.samples$count.tissue[which(count.samples$class=="nepc" & count.samples$count.tissue>1)]),sep="|")))

text(1.5,0.86,paste("n.points = ",paste(length(dist.adeno),length(dist.nepc),sep="|")))


boxplot(dist.adeno.plasma_tissue,dist.nepc.plasma_tissue,varwidth = T,names = c("CRPC-Adeno","CRPC-NE"),main="Fraction of plasma in tissue",ylab="SNV Similarity (Non-synonymous)",ylim=c(0,1),outline = F) #ylim=c(0,0.3)
stripchart(list(dist.adeno.plasma_tissue,dist.nepc.plasma_tissue), vertical=T, method="jitter",add=TRUE,pch=19,col=c(crpc.color.rgb,nepc.color.rgb),jitter = 0.2)
p = wilcox.test(x = dist.adeno.plasma_tissue,dist.nepc.plasma_tissue)
p = formatC(p$p.value,format="e",digits=2)
text(1.5,1,paste("p = ",p,sep=""))

text(1.5,1*0.96,paste("n.patients =",paste(length(unique(count.samples$id.patient[which(count.samples$class=="adeno")])),
                                                  length(unique(count.samples$id.patient[which(count.samples$class=="nepc")])),sep="|")))

text(1.5,1*0.91,paste("n.plasma =",paste(sum(count.samples$count.plasma[which(count.samples$class=="adeno")]),
                                                sum(count.samples$count.plasma[which(count.samples$class=="nepc")]),sep="|")))

text(1.5,1*0.86,paste("n.tissues =",paste(sum(count.samples$count.tissue[which(count.samples$class=="adeno")]),
                                                 sum(count.samples$count.tissue[which(count.samples$class=="nepc")]),sep="|")))

text(1.5,1*0.81,paste("n.points = ",paste(length(dist.adeno.plasma_tissue),length(dist.nepc.plasma_tissue),sep="|")))

boxplot(dist.adeno.tissue_plasma,dist.nepc.tissue_plasma,varwidth = T,names = c("CRPC-Adeno","CRPC-NE"),main="Fraction of tissue in plasma",ylab="SNV Similarity (Non-synonymous)",ylim=c(0,1),outline = F)
stripchart(list(dist.adeno.tissue_plasma,dist.nepc.tissue_plasma), vertical=T, method="jitter",add=TRUE,pch=19,col=c(crpc.color.rgb,nepc.color.rgb),jitter = 0.2)
p = wilcox.test(x = dist.adeno.tissue_plasma,dist.nepc.tissue_plasma)
p = formatC(p$p.value,format="e",digits=2)
text(1.5,1,paste("p = ",p,sep=""))

text(1.5,0.96,paste("n.patients =",paste(length(unique(count.samples$id.patient[which(count.samples$class=="adeno")])),
                                         length(unique(count.samples$id.patient[which(count.samples$class=="nepc")])),sep="|")))

text(1.5,0.91,paste("n.plasma =",paste(sum(count.samples$count.plasma[which(count.samples$class=="adeno")]),
                                       sum(count.samples$count.plasma[which(count.samples$class=="nepc")]),sep="|")))

text(1.5,0.86,paste("n.tissues =",paste(sum(count.samples$count.tissue[which(count.samples$class=="adeno")]),
                                        sum(count.samples$count.tissue[which(count.samples$class=="nepc")]),sep="|")))

text(1.5,0.81,paste("n.points = ",paste(length(dist.adeno.tissue_plasma),length(dist.nepc.tissue_plasma),sep="|")))

dev.off()
