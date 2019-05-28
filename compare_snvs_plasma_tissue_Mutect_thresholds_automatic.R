library(reshape2)
library(ggplot2)
library(stringr)

min.tc = 0

min.af.case = 0.05
min.cov.case = 10

max.af.ctrl = 0.01
min.cov.ctrl = 0

setwd("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/")

# Plasma sif 
plasma.sif = read.delim("/elaborazioni/sharedCO/Abemus_data_analysis/Samples_info_files/IPM_PLASMA/samples_70_plasma.tsv",as.is = T,header = T,stringsAsFactors = F)
plasma.sif$id.patient = gsub(plasma.sif$plasma,pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = "")
plasma.sif$type = "plasma"

# Tissue sif HALO
halo.sif = read.delim("/elaborazioni/sharedCO/Abemus_data_analysis/Samples_info_files/IPM/sif_tissue_samples_haloplex.tsv",as.is = T,header = T,stringsAsFactors = F)
ids = c()
for(i in 1:nrow(halo.sif)){
  ids = c(ids,unlist(strsplit(halo.sif$patient[i],split = "_"))[1])
}
halo.sif$id.patient = ids
halo.sif$type = "halo"

# Tissue sif SureSelect
sure.sif = read.delim("/elaborazioni/sharedCO/Abemus_data_analysis/Samples_info_files/IPM/sif_tissue_samples_sureselect_COMPLETE.tsv",as.is = T,header = T,stringsAsFactors = F)
ids = c()
for(i in 1:nrow(sure.sif)){
  ids = c(ids,unlist(strsplit(sure.sif$patient[i],split = "_"))[1])
}
sure.sif$id.patient = ids
sure.sif$type = "sure"

### Mutect version 
# Mutect calls plasma
calls = readLines(con = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect/Plasma/calls_mutect.txt")

plasma.calls = data.frame(id.case = str_split_fixed(calls, "/", 12)[,11], calls = calls, stringsAsFactors = F)
plasma.calls$id.patient = gsub(plasma.calls$id.case, pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = "")

# Strelka2 calls halo
halo.calls = readLines(con = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect/Tissue_HaloPlex/calls_mutect.txt")

halo.calls = data.frame(id.patient = str_split_fixed(halo.calls, "/", 12)[,11], calls = halo.calls, stringsAsFactors = F)

# Strelka2 calls sure
sure.calls = readLines(con = "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect/Tissue_SureSelect/calls_mutect.txt")

sure.calls = data.frame(id.patient = str_split_fixed(sure.calls, "/", 12)[,11], calls = sure.calls, stringsAsFactors = F)

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

# Genelist
goi = read.delim("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/abemus_plasma_tissues_distances/Genes_List_Plasma_Study_180321.txt",as.is = T,stringsAsFactors = F)

# combinations
patients = unique(gsub(anns.plasma,pattern = "-1st|-TP1|-TP2|-TP3|-TP4",replacement = ""))
length(patients)



for (min.cov.case in c(10,20,30,50)) {
  if(min.tc != 0){
    dir <- file.path("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/", paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case, "_mintc_", min.tc))
    dir.create(path = dir)
    setwd(dir)
  }else{
    dir <- file.path("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/", paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case))
    dir.create(path = dir)
    setwd(dir)
  }
  
  pdf("plots_and_distances_matrices.pdf")
  
  for(p in patients){
    message(p)
    x = c(plasma.calls$calls[which(plasma.calls$id.patient==p)],
          halo.calls$calls[which(halo.calls$id.patient==p)],
          sure.calls$calls[which(sure.calls$id.patient==p)])
    basename(x)
    if(length(x)==1){
      message(p," only plasma")
      next
    } else {
      pairs=combn(x,2,simplify = T)
      pairs.name = pairs
      min.purity = c()
      for(j in 1:ncol(pairs.name)){
        pairs.name[1,j] = str_split_fixed(pairs.name[1,j], "/", n = 12)[,11] #9 instead of 12 if in local
        pairs.name[2,j] = str_split_fixed(pairs.name[2,j], "/", n = 12)[,11]
        pairs.name[1,j] = gsub(pairs.name[1,j], pattern = "_HALO|-HALO", replacement = "")
        pairs.name[2,j] = gsub(pairs.name[2,j], pattern = "_HALO|-HALO", replacement = "")
        pairs.name[which(pairs.name==p,arr.ind = T)] = paste0(pairs.name[which(pairs.name==p,arr.ind = T)],"-1st")
        min.purity = c(min.purity,min(gadm$purity[which(gadm$sample==pairs.name[1,j])],gadm$purity[which(gadm$sample==pairs.name[2,j])]))
      }
      # filter based on tumor content
      pairs = pairs[,which(min.purity >= min.tc),drop=F]
      pairs.name = pairs.name[,which(min.purity >= min.tc),drop=F]
      
      if(ncol(pairs)==0){
        message("no samples with enough tumor purity")
        next
      }
      
      dimnames = unique(c(pairs.name[1,],pairs.name[2,]))
      dimnames = data.frame(sample=dimnames,stringsAsFactors = F)
      dimnames = merge(x = dimnames,y = anns,by = "sample",all.x = T)
      dimnames$date = as.Date(dimnames$date,format = "%m/%d/%Y")
      dimnames = dimnames[with(dimnames, order(date,decreasing = F)), ]
      
      plasma.1st = paste0(p,"-1st") 
      keep.plasma.first = c(which(dimnames$sample==plasma.1st),which(dimnames$sample!=plasma.1st))
      dimnames = dimnames$sample[keep.plasma.first]
      dimnames = gsub(dimnames,pattern = "-1st",replacement = "")
      
      m = matrix(data = 1,nrow = length(dimnames),ncol = length(dimnames),dimnames = list(dimnames,dimnames))
      dd = c()
      for(i in 1:ncol(pairs)){
        # sample a
        a=read.delim(file = pairs[1,i],header = T,stringsAsFactors = F)
        a <- a[which(a$af_case > min.af.case & a$cov_case > min.cov.case),]
        a <- a[which(a$af_control < max.af.ctrl & a$cov_control > min.cov.ctrl),]
        # at least 2 alternatives
        a <- a[which(a$af_case*a$cov_case >= 2),]
        # sample b
        b=read.delim(file = pairs[2,i],header = T,stringsAsFactors = F)
        b <- b[which(b$af_case > min.af.case & b$cov_case > min.cov.case),]
        b <- b[which(b$af_control < max.af.ctrl & b$cov_control > min.cov.ctrl),]
        b <- b[which(b$af_case*b$cov_case >= 2),]
        
        shared = merge(x = a,y = b,by = "group",all = F,suffixes = c("_a","_b"))
        #p1.name = clean.name(x = pairs[1,i],suffixes = suffixes)
        p1.name = str_split_fixed(pairs[1,i], "/", n = 12)[,11]
        p1.name = gsub(p1.name, pattern = "_HALO|-HALO", replacement = "")
        #p2.name = clean.name(x = pairs[2,i],suffixes = suffixes)
        p2.name = str_split_fixed(pairs[2,i], "/", n = 12)[,11]
        p2.name = gsub(p2.name, pattern = "_HALO|-HALO", replacement = "")
        #m[p1.name,p2.name] <- length(unique(shared$group))/length(unique(a$group))
        #m[p2.name,p1.name] <- length(unique(shared$group))/length(unique(b$group))
        m[p2.name,p1.name] <- length(unique(shared$group))/length(unique(a$group))
        m[p1.name,p2.name] <- length(unique(shared$group))/length(unique(b$group))
        
        td1 = data.frame(a=p1.name,b=p2.name,similarity=m[p2.name,p1.name],stringsAsFactors = F)
        td2 = data.frame(a=p2.name,b=p1.name,similarity=m[p1.name,p2.name],stringsAsFactors = F)
        dd = rbind(dd,td1,td2)
        
        ids.shared = unique(shared$group)
        ids.only.a = setdiff(unique(a$group),ids.shared)
        ids.only.b = setdiff(unique(b$group),ids.shared)
        
        only.a = a[which(a$group%in%ids.only.a),]
        only.b = b[which(b$group%in%ids.only.b),]
        
        par(pty="s")
        
        col.shared = rgb(red = 49,green = 130,blue = 189,maxColorValue = 255,alpha = 150)
        col.a = rgb(217,95,14,maxColorValue = 255,alpha = 150)
        col.b = rgb(49,163,84,maxColorValue = 255,alpha = 150)
        
        plot(x = shared$af_case_b,y = shared$af_case_a,xlim = c(-0.05,1.0),ylim = c(-0.05,1.0),pch=16,main=p,xlab=paste0("Uncorrected AF in ",p2.name),ylab=paste0("Uncorrected AF in ",p1.name),col=col.shared)
        grid()
        abline(h = 0,v = 0,col="grey20")
        abline(coef = c(0,1),col="grey80",lwd=1,lty="dashed")
        x.coords = rnorm(nrow(only.a),mean = -0.05,sd = 0.007)
        lines(y = only.a$af_case,x = x.coords,ylim=c(0,1),xlim=c(0,1),type = "p",pch=16,col=col.a)
        y.coords = rnorm(nrow(only.b),mean = -0.05,sd = 0.007)
        lines(y = y.coords,x = only.b$af_case,ylim=c(0,1),xlim=c(0,1),type = "p",pch=16,col=col.b)
        
        # all Hugo names 
        if(F & length(shared$Hugo_Symbol_a != 0)){
          #text(y = only.a$af_case,x = x.coords,labels = only.a$Hugo_Symbol,pos = 3,cex = 0.5)
          if(nrow(shared)>0){text(y = shared$af_case_a,x = shared$af_case_b,labels = shared$Hugo_Symbol_a,pos = 3,cex = 0.5)}
          #text(y = y.coords,x = only.b$af_case,labels = only.b$Hugo_Symbol,pos = 3,cex = 0.5)
        }
        # only Hugo names of interest
        if(T){
          shared$Hugo_Symbol_a[which(!shared$Hugo_Symbol_a%in%goi$symbol)] <- NA
          only.a$Hugo_Symbol[which(!only.a$Hugo_Symbol%in%goi$symbol)] <- NA
          only.b$Hugo_Symbol[which(!only.b$Hugo_Symbol%in%goi$symbol)] <- NA
          if(nrow(shared)>0){text(y = shared$af_case_a,x = shared$af_case_b,labels = shared$Hugo_Symbol_a,pos = 3,cex = 0.5)}
          if(nrow(only.a)>0){text(y = only.a$af_case,x = x.coords,labels = only.a$Hugo_Symbol,pos = 3,cex = 0.5)}
          if(nrow(only.b)>0){text(y = y.coords,x = only.b$af_case,labels = only.b$Hugo_Symbol,pos = 3,cex = 0.5)}
        }
        legend(box.lty=0,x = 0.05,y = 1,legend=c(paste0("n = ",nrow(shared)),paste0("n = ",nrow(only.a)),paste0("n = ",nrow(only.b))),pch = 16,col=c(col.shared,col.a,col.b), cex=0.9)
      }
      save(m,dd,file = paste0("distances_matrix_",p,".RData"),compress = T)
    }
    dat = m
    dat2 = melt(dat)
    pp = ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
      geom_tile(aes(fill = value)) + 
      geom_text(aes(label = round(dat2$value, 3))) +
      scale_fill_gradientn(limits = c(0,1),colours = c("#e5f5e0","#a1d99b","#31a354")) +
      ylab("") + xlab("") + coord_fixed() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)) +
      ggtitle(paste("Fraction of included SNVs ",p,sep=""))
    
    print(pp)
    
  }
  
  dev.off()
  
  
  ### all labels
  pdf("plots_and_distances_matrices_sharedlabels.pdf")
  
  for(p in patients){
    message(p)
    x = c(plasma.calls$calls[which(plasma.calls$id.patient==p)],
          halo.calls$calls[which(halo.calls$id.patient==p)],
          sure.calls$calls[which(sure.calls$id.patient==p)])
    basename(x)
    if(length(x)==1){
      message(p," only plasma")
      next
    } else {
      pairs=combn(x,2,simplify = T)
      pairs.name = pairs
      min.purity = c()
      for(j in 1:ncol(pairs.name)){
        pairs.name[1,j] = str_split_fixed(pairs.name[1,j], "/", n = 12)[,11] #9 instead of 12 if in local
        pairs.name[2,j] = str_split_fixed(pairs.name[2,j], "/", n = 12)[,11]
        pairs.name[1,j] = gsub(pairs.name[1,j], pattern = "_HALO|-HALO", replacement = "")
        pairs.name[2,j] = gsub(pairs.name[2,j], pattern = "_HALO|-HALO", replacement = "")
        pairs.name[which(pairs.name==p,arr.ind = T)] = paste0(pairs.name[which(pairs.name==p,arr.ind = T)],"-1st")
        min.purity = c(min.purity,min(gadm$purity[which(gadm$sample==pairs.name[1,j])],gadm$purity[which(gadm$sample==pairs.name[2,j])]))
      }
      # filter based on tumor content
      pairs = pairs[,which(min.purity >= min.tc),drop=F]
      pairs.name = pairs.name[,which(min.purity >= min.tc),drop=F]
      
      if(ncol(pairs)==0){
        message("no samples with enough tumor purity")
        next
      }
      
      dimnames = unique(c(pairs.name[1,],pairs.name[2,]))
      dimnames = data.frame(sample=dimnames,stringsAsFactors = F)
      dimnames = merge(x = dimnames,y = anns,by = "sample",all.x = T)
      dimnames$date = as.Date(dimnames$date,format = "%m/%d/%Y")
      dimnames = dimnames[with(dimnames, order(date,decreasing = F)), ]
      
      plasma.1st = paste0(p,"-1st") 
      keep.plasma.first = c(which(dimnames$sample==plasma.1st),which(dimnames$sample!=plasma.1st))
      dimnames = dimnames$sample[keep.plasma.first]
      dimnames = gsub(dimnames,pattern = "-1st",replacement = "")
      
      m = matrix(data = 1,nrow = length(dimnames),ncol = length(dimnames),dimnames = list(dimnames,dimnames))
      dd = c()
      for(i in 1:ncol(pairs)){
        # sample a
        a=read.delim(file = pairs[1,i],header = T,stringsAsFactors = F)
        a <- a[which(a$af_case > min.af.case & a$cov_case > min.cov.case),]
        a <- a[which(a$af_control < max.af.ctrl & a$cov_control > min.cov.ctrl),]
        # at least 2 alternatives
        a <- a[which(a$af_case*a$cov_case >= 2),]
        # sample b
        b=read.delim(file = pairs[2,i],header = T,stringsAsFactors = F)
        b <- b[which(b$af_case > min.af.case & b$cov_case > min.cov.case),]
        b <- b[which(b$af_control < max.af.ctrl & b$cov_control > min.cov.ctrl),]
        b <- b[which(b$af_case*b$cov_case >= 2),]
        
        shared = merge(x = a,y = b,by = "group",all = F,suffixes = c("_a","_b"))
        #p1.name = clean.name(x = pairs[1,i],suffixes = suffixes)
        p1.name = str_split_fixed(pairs[1,i], "/", n = 12)[,11]
        p1.name = gsub(p1.name, pattern = "_HALO|-HALO", replacement = "")
        #p2.name = clean.name(x = pairs[2,i],suffixes = suffixes)
        p2.name = str_split_fixed(pairs[2,i], "/", n = 12)[,11]
        p2.name = gsub(p2.name, pattern = "_HALO|-HALO", replacement = "")
        
        #m[p1.name,p2.name] <- length(unique(shared$group))/length(unique(a$group))
        #m[p2.name,p1.name] <- length(unique(shared$group))/length(unique(b$group))
        m[p2.name,p1.name] <- length(unique(shared$group))/length(unique(a$group))
        m[p1.name,p2.name] <- length(unique(shared$group))/length(unique(b$group))
        
        td1 = data.frame(a=p1.name,b=p2.name,similarity=m[p2.name,p1.name],stringsAsFactors = F)
        td2 = data.frame(a=p2.name,b=p1.name,similarity=m[p1.name,p2.name],stringsAsFactors = F)
        dd = rbind(dd,td1,td2)
        
        ids.shared = unique(shared$group)
        ids.only.a = setdiff(unique(a$group),ids.shared)
        ids.only.b = setdiff(unique(b$group),ids.shared)
        
        only.a = a[which(a$group%in%ids.only.a),]
        only.b = b[which(b$group%in%ids.only.b),]
        
        par(pty="s")
        
        col.shared = rgb(red = 49,green = 130,blue = 189,maxColorValue = 255,alpha = 150)
        col.a = rgb(217,95,14,maxColorValue = 255,alpha = 150)
        col.b = rgb(49,163,84,maxColorValue = 255,alpha = 150)
        
        plot(x = shared$af_case_b,y = shared$af_case_a,xlim = c(-0.05,1.0),ylim = c(-0.05,1.0),pch=16,main=p,xlab=paste0("Uncorrected AF in ",p2.name),ylab=paste0("Uncorrected AF in ",p1.name),col=col.shared)
        grid()
        abline(h = 0,v = 0,col="grey20")
        abline(coef = c(0,1),col="grey80",lwd=1,lty="dashed")
        x.coords = rnorm(nrow(only.a),mean = -0.05,sd = 0.007)
        lines(y = only.a$af_case,x = x.coords,ylim=c(0,1),xlim=c(0,1),type = "p",pch=16,col=col.a)
        y.coords = rnorm(nrow(only.b),mean = -0.05,sd = 0.007)
        lines(y = y.coords,x = only.b$af_case,ylim=c(0,1),xlim=c(0,1),type = "p",pch=16,col=col.b)
        
        # all Hugo names 
        if(T & length(shared$Hugo_Symbol_a != 0)){
          #text(y = only.a$af_case,x = x.coords,labels = only.a$Hugo_Symbol,pos = 3,cex = 0.5)
          if(nrow(shared)>0){text(y = shared$af_case_a,x = shared$af_case_b,labels = shared$Hugo_Symbol_a,pos = 3,cex = 0.5)}
          #text(y = y.coords,x = only.b$af_case,labels = only.b$Hugo_Symbol,pos = 3,cex = 0.5)
        }
        # only Hugo names of interest
        if(F){
          shared$Hugo_Symbol_a[which(!shared$Hugo_Symbol_a%in%goi$symbol)] <- NA
          only.a$Hugo_Symbol[which(!only.a$Hugo_Symbol%in%goi$symbol)] <- NA
          only.b$Hugo_Symbol[which(!only.b$Hugo_Symbol%in%goi$symbol)] <- NA
          if(nrow(shared)>0){text(y = shared$af_case_a,x = shared$af_case_b,labels = shared$Hugo_Symbol_a,pos = 3,cex = 0.5)}
          if(nrow(only.a)>0){text(y = only.a$af_case,x = x.coords,labels = only.a$Hugo_Symbol,pos = 3,cex = 0.5)}
          if(nrow(only.b)>0){text(y = y.coords,x = only.b$af_case,labels = only.b$Hugo_Symbol,pos = 3,cex = 0.5)}
        }
        legend(box.lty=0,x = 0.05,y = 1,legend=c(paste0("n = ",nrow(shared)),paste0("n = ",nrow(only.a)),paste0("n = ",nrow(only.b))),pch = 16,col=c(col.shared,col.a,col.b), cex=0.9)
      }
      #save(m,dd,file = paste0("distances_matrix_",p,".RData"),compress = T)
    }
    dat = m
    dat2 = melt(dat)
    pp = ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
      geom_tile(aes(fill = value)) + 
      geom_text(aes(label = round(dat2$value, 3))) +
      scale_fill_gradientn(limits = c(0,1),colours = c("#e5f5e0","#a1d99b","#31a354")) +
      ylab("") + xlab("") + coord_fixed() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)) +
      ggtitle(paste("Fraction of included SNVs ",p,sep=""))
    
    print(pp)
    
  }
  
  dev.off()
}


min.af.case = 0.01
min.cov.case = 10

if(min.tc != 0){
  dir <- file.path("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/", paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case, "_mintc_", min.tc))
  dir.create(path = dir)
  setwd(dir)
}else{
  dir <- file.path("/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_plasma_tissue_distances/", paste0("afn_",max.af.ctrl, "_aft_", min.af.case, "_mincov_", min.cov.case))
  dir.create(path = dir)
  setwd(dir)
}

pdf("plots_and_distances_matrices.pdf")

for(p in patients){
  message(p)
  x = c(plasma.calls$calls[which(plasma.calls$id.patient==p)],
        halo.calls$calls[which(halo.calls$id.patient==p)],
        sure.calls$calls[which(sure.calls$id.patient==p)])
  basename(x)
  if(length(x)==1){
    message(p," only plasma")
    next
  } else {
    pairs=combn(x,2,simplify = T)
    pairs.name = pairs
    min.purity = c()
    for(j in 1:ncol(pairs.name)){
      pairs.name[1,j] = str_split_fixed(pairs.name[1,j], "/", n = 12)[,11] #9 instead of 12 if in local
      pairs.name[2,j] = str_split_fixed(pairs.name[2,j], "/", n = 12)[,11]
      pairs.name[1,j] = gsub(pairs.name[1,j], pattern = "_HALO|-HALO", replacement = "")
      pairs.name[2,j] = gsub(pairs.name[2,j], pattern = "_HALO|-HALO", replacement = "")
      pairs.name[which(pairs.name==p,arr.ind = T)] = paste0(pairs.name[which(pairs.name==p,arr.ind = T)],"-1st")
      min.purity = c(min.purity,min(gadm$purity[which(gadm$sample==pairs.name[1,j])],gadm$purity[which(gadm$sample==pairs.name[2,j])]))
    }
    # filter based on tumor content
    pairs = pairs[,which(min.purity >= min.tc),drop=F]
    pairs.name = pairs.name[,which(min.purity >= min.tc),drop=F]
    
    if(ncol(pairs)==0){
      message("no samples with enough tumor purity")
      next
    }
    
    dimnames = unique(c(pairs.name[1,],pairs.name[2,]))
    dimnames = data.frame(sample=dimnames,stringsAsFactors = F)
    dimnames = merge(x = dimnames,y = anns,by = "sample",all.x = T)
    dimnames$date = as.Date(dimnames$date,format = "%m/%d/%Y")
    dimnames = dimnames[with(dimnames, order(date,decreasing = F)), ]
    
    plasma.1st = paste0(p,"-1st") 
    keep.plasma.first = c(which(dimnames$sample==plasma.1st),which(dimnames$sample!=plasma.1st))
    dimnames = dimnames$sample[keep.plasma.first]
    dimnames = gsub(dimnames,pattern = "-1st",replacement = "")
    
    m = matrix(data = 1,nrow = length(dimnames),ncol = length(dimnames),dimnames = list(dimnames,dimnames))
    dd = c()
    for(i in 1:ncol(pairs)){
      # sample a
      a=read.delim(file = pairs[1,i],header = T,stringsAsFactors = F)
      a <- a[which(a$af_case > min.af.case & a$cov_case > min.cov.case),]
      a <- a[which(a$af_control < max.af.ctrl & a$cov_control > min.cov.ctrl),]
      # at least 2 alternatives
      a <- a[which(a$af_case*a$cov_case >= 2),]
      # sample b
      b=read.delim(file = pairs[2,i],header = T,stringsAsFactors = F)
      b <- b[which(b$af_case > min.af.case & b$cov_case > min.cov.case),]
      b <- b[which(b$af_control < max.af.ctrl & b$cov_control > min.cov.ctrl),]
      b <- b[which(b$af_case*b$cov_case >= 2),]
      
      shared = merge(x = a,y = b,by = "group",all = F,suffixes = c("_a","_b"))
      #p1.name = clean.name(x = pairs[1,i],suffixes = suffixes)
      p1.name = str_split_fixed(pairs[1,i], "/", n = 12)[,11]
      p1.name = gsub(p1.name, pattern = "_HALO|-HALO", replacement = "")
      #p2.name = clean.name(x = pairs[2,i],suffixes = suffixes)
      p2.name = str_split_fixed(pairs[2,i], "/", n = 12)[,11]
      p2.name = gsub(p2.name, pattern = "_HALO|-HALO", replacement = "")
      
      #m[p1.name,p2.name] <- length(unique(shared$group))/length(unique(a$group))
      #m[p2.name,p1.name] <- length(unique(shared$group))/length(unique(b$group))
      m[p2.name,p1.name] <- length(unique(shared$group))/length(unique(a$group))
      m[p1.name,p2.name] <- length(unique(shared$group))/length(unique(b$group))
      
      td1 = data.frame(a=p1.name,b=p2.name,similarity=m[p2.name,p1.name],stringsAsFactors = F)
      td2 = data.frame(a=p2.name,b=p1.name,similarity=m[p1.name,p2.name],stringsAsFactors = F)
      dd = rbind(dd,td1,td2)
      
      ids.shared = unique(shared$group)
      ids.only.a = setdiff(unique(a$group),ids.shared)
      ids.only.b = setdiff(unique(b$group),ids.shared)
      
      only.a = a[which(a$group%in%ids.only.a),]
      only.b = b[which(b$group%in%ids.only.b),]
      
      par(pty="s")
      
      col.shared = rgb(red = 49,green = 130,blue = 189,maxColorValue = 255,alpha = 150)
      col.a = rgb(217,95,14,maxColorValue = 255,alpha = 150)
      col.b = rgb(49,163,84,maxColorValue = 255,alpha = 150)
      
      plot(x = shared$af_case_b,y = shared$af_case_a,xlim = c(-0.05,1.0),ylim = c(-0.05,1.0),pch=16,main=p,xlab=paste0("Uncorrected AF in ",p2.name),ylab=paste0("Uncorrected AF in ",p1.name),col=col.shared)
      grid()
      abline(h = 0,v = 0,col="grey20")
      abline(coef = c(0,1),col="grey80",lwd=1,lty="dashed")
      x.coords = rnorm(nrow(only.a),mean = -0.05,sd = 0.007)
      lines(y = only.a$af_case,x = x.coords,ylim=c(0,1),xlim=c(0,1),type = "p",pch=16,col=col.a)
      y.coords = rnorm(nrow(only.b),mean = -0.05,sd = 0.007)
      lines(y = y.coords,x = only.b$af_case,ylim=c(0,1),xlim=c(0,1),type = "p",pch=16,col=col.b)
      
      # all Hugo names 
      if(F & length(shared$Hugo_Symbol_a != 0)){
        #text(y = only.a$af_case,x = x.coords,labels = only.a$Hugo_Symbol,pos = 3,cex = 0.5)
        if(nrow(shared)>0){text(y = shared$af_case_a,x = shared$af_case_b,labels = shared$Hugo_Symbol_a,pos = 3,cex = 0.5)}
        #text(y = y.coords,x = only.b$af_case,labels = only.b$Hugo_Symbol,pos = 3,cex = 0.5)
      }
      # only Hugo names of interest
      if(T){
        shared$Hugo_Symbol_a[which(!shared$Hugo_Symbol_a%in%goi$symbol)] <- NA
        only.a$Hugo_Symbol[which(!only.a$Hugo_Symbol%in%goi$symbol)] <- NA
        only.b$Hugo_Symbol[which(!only.b$Hugo_Symbol%in%goi$symbol)] <- NA
        if(nrow(shared)>0){text(y = shared$af_case_a,x = shared$af_case_b,labels = shared$Hugo_Symbol_a,pos = 3,cex = 0.5)}
        if(nrow(only.a)>0){text(y = only.a$af_case,x = x.coords,labels = only.a$Hugo_Symbol,pos = 3,cex = 0.5)}
        if(nrow(only.b)>0){text(y = y.coords,x = only.b$af_case,labels = only.b$Hugo_Symbol,pos = 3,cex = 0.5)}
      }
      legend(box.lty=0,x = 0.05,y = 1,legend=c(paste0("n = ",nrow(shared)),paste0("n = ",nrow(only.a)),paste0("n = ",nrow(only.b))),pch = 16,col=c(col.shared,col.a,col.b), cex=0.9)
    }
    save(m,dd,file = paste0("distances_matrix_",p,".RData"),compress = T)
  }
  dat = m
  dat2 = melt(dat)
  pp = ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
    geom_tile(aes(fill = value)) + 
    geom_text(aes(label = round(dat2$value, 3))) +
    scale_fill_gradientn(limits = c(0,1),colours = c("#e5f5e0","#a1d99b","#31a354")) +
    ylab("") + xlab("") + coord_fixed() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)) +
    ggtitle(paste("Fraction of included SNVs ",p,sep=""))
  
  print(pp)
  
}

dev.off()


### all labels
pdf("plots_and_distances_matrices_sharedlabels.pdf")

for(p in patients){
  message(p)
  x = c(plasma.calls$calls[which(plasma.calls$id.patient==p)],
        halo.calls$calls[which(halo.calls$id.patient==p)],
        sure.calls$calls[which(sure.calls$id.patient==p)])
  basename(x)
  if(length(x)==1){
    message(p," only plasma")
    next
  } else {
    pairs=combn(x,2,simplify = T)
    pairs.name = pairs
    min.purity = c()
    for(j in 1:ncol(pairs.name)){
      pairs.name[1,j] = str_split_fixed(pairs.name[1,j], "/", n = 12)[,11] #9 instead of 12 if in local
      pairs.name[2,j] = str_split_fixed(pairs.name[2,j], "/", n = 12)[,11]
      pairs.name[1,j] = gsub(pairs.name[1,j], pattern = "_HALO|-HALO", replacement = "")
      pairs.name[2,j] = gsub(pairs.name[2,j], pattern = "_HALO|-HALO", replacement = "")
      pairs.name[which(pairs.name==p,arr.ind = T)] = paste0(pairs.name[which(pairs.name==p,arr.ind = T)],"-1st")
      min.purity = c(min.purity,min(gadm$purity[which(gadm$sample==pairs.name[1,j])],gadm$purity[which(gadm$sample==pairs.name[2,j])]))
    }
    # filter based on tumor content
    pairs = pairs[,which(min.purity >= min.tc),drop=F]
    pairs.name = pairs.name[,which(min.purity >= min.tc),drop=F]
    
    if(ncol(pairs)==0){
      message("no samples with enough tumor purity")
      next
    }
    
    dimnames = unique(c(pairs.name[1,],pairs.name[2,]))
    dimnames = data.frame(sample=dimnames,stringsAsFactors = F)
    dimnames = merge(x = dimnames,y = anns,by = "sample",all.x = T)
    dimnames$date = as.Date(dimnames$date,format = "%m/%d/%Y")
    dimnames = dimnames[with(dimnames, order(date,decreasing = F)), ]
    
    plasma.1st = paste0(p,"-1st") 
    keep.plasma.first = c(which(dimnames$sample==plasma.1st),which(dimnames$sample!=plasma.1st))
    dimnames = dimnames$sample[keep.plasma.first]
    dimnames = gsub(dimnames,pattern = "-1st",replacement = "")
    
    m = matrix(data = 1,nrow = length(dimnames),ncol = length(dimnames),dimnames = list(dimnames,dimnames))
    dd = c()
    for(i in 1:ncol(pairs)){
      # sample a
      a=read.delim(file = pairs[1,i],header = T,stringsAsFactors = F)
      a <- a[which(a$af_case > min.af.case & a$cov_case > min.cov.case),]
      a <- a[which(a$af_control < max.af.ctrl & a$cov_control > min.cov.ctrl),]
      # at least 2 alternatives
      a <- a[which(a$af_case*a$cov_case >= 2),]
      # sample b
      b=read.delim(file = pairs[2,i],header = T,stringsAsFactors = F)
      b <- b[which(b$af_case > min.af.case & b$cov_case > min.cov.case),]
      b <- b[which(b$af_control < max.af.ctrl & b$cov_control > min.cov.ctrl),]
      b <- b[which(b$af_case*b$cov_case >= 2),]
      
      shared = merge(x = a,y = b,by = "group",all = F,suffixes = c("_a","_b"))
      #p1.name = clean.name(x = pairs[1,i],suffixes = suffixes)
      p1.name = str_split_fixed(pairs[1,i], "/", n = 12)[,11]
      p1.name = gsub(p1.name, pattern = "_HALO|-HALO", replacement = "")
      #p2.name = clean.name(x = pairs[2,i],suffixes = suffixes)
      p2.name = str_split_fixed(pairs[2,i], "/", n = 12)[,11]
      p2.name = gsub(p2.name, pattern = "_HALO|-HALO", replacement = "")
      
      #m[p1.name,p2.name] <- length(unique(shared$group))/length(unique(a$group))
      #m[p2.name,p1.name] <- length(unique(shared$group))/length(unique(b$group))
      m[p2.name,p1.name] <- length(unique(shared$group))/length(unique(a$group))
      m[p1.name,p2.name] <- length(unique(shared$group))/length(unique(b$group))
      
      td1 = data.frame(a=p1.name,b=p2.name,similarity=m[p2.name,p1.name],stringsAsFactors = F)
      td2 = data.frame(a=p2.name,b=p1.name,similarity=m[p1.name,p2.name],stringsAsFactors = F)
      dd = rbind(dd,td1,td2)
      
      ids.shared = unique(shared$group)
      ids.only.a = setdiff(unique(a$group),ids.shared)
      ids.only.b = setdiff(unique(b$group),ids.shared)
      
      only.a = a[which(a$group%in%ids.only.a),]
      only.b = b[which(b$group%in%ids.only.b),]
      
      par(pty="s")
      
      col.shared = rgb(red = 49,green = 130,blue = 189,maxColorValue = 255,alpha = 150)
      col.a = rgb(217,95,14,maxColorValue = 255,alpha = 150)
      col.b = rgb(49,163,84,maxColorValue = 255,alpha = 150)
      
      plot(x = shared$af_case_b,y = shared$af_case_a,xlim = c(-0.05,1.0),ylim = c(-0.05,1.0),pch=16,main=p,xlab=paste0("Uncorrected AF in ",p2.name),ylab=paste0("Uncorrected AF in ",p1.name),col=col.shared)
      grid()
      abline(h = 0,v = 0,col="grey20")
      abline(coef = c(0,1),col="grey80",lwd=1,lty="dashed")
      x.coords = rnorm(nrow(only.a),mean = -0.05,sd = 0.007)
      lines(y = only.a$af_case,x = x.coords,ylim=c(0,1),xlim=c(0,1),type = "p",pch=16,col=col.a)
      y.coords = rnorm(nrow(only.b),mean = -0.05,sd = 0.007)
      lines(y = y.coords,x = only.b$af_case,ylim=c(0,1),xlim=c(0,1),type = "p",pch=16,col=col.b)
      
      # all Hugo names 
      if(T & length(shared$Hugo_Symbol_a != 0)){
        #text(y = only.a$af_case,x = x.coords,labels = only.a$Hugo_Symbol,pos = 3,cex = 0.5)
        if(nrow(shared)>0){text(y = shared$af_case_a,x = shared$af_case_b,labels = shared$Hugo_Symbol_a,pos = 3,cex = 0.5)}
        #text(y = y.coords,x = only.b$af_case,labels = only.b$Hugo_Symbol,pos = 3,cex = 0.5)
      }
      # only Hugo names of interest
      if(F){
        shared$Hugo_Symbol_a[which(!shared$Hugo_Symbol_a%in%goi$symbol)] <- NA
        only.a$Hugo_Symbol[which(!only.a$Hugo_Symbol%in%goi$symbol)] <- NA
        only.b$Hugo_Symbol[which(!only.b$Hugo_Symbol%in%goi$symbol)] <- NA
        if(nrow(shared)>0){text(y = shared$af_case_a,x = shared$af_case_b,labels = shared$Hugo_Symbol_a,pos = 3,cex = 0.5)}
        if(nrow(only.a)>0){text(y = only.a$af_case,x = x.coords,labels = only.a$Hugo_Symbol,pos = 3,cex = 0.5)}
        if(nrow(only.b)>0){text(y = y.coords,x = only.b$af_case,labels = only.b$Hugo_Symbol,pos = 3,cex = 0.5)}
      }
      legend(box.lty=0,x = 0.05,y = 1,legend=c(paste0("n = ",nrow(shared)),paste0("n = ",nrow(only.a)),paste0("n = ",nrow(only.b))),pch = 16,col=c(col.shared,col.a,col.b), cex=0.9)
    }
    #save(m,dd,file = paste0("distances_matrix_",p,".RData"),compress = T)
  }
  dat = m
  dat2 = melt(dat)
  pp = ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
    geom_tile(aes(fill = value)) + 
    geom_text(aes(label = round(dat2$value, 3))) +
    scale_fill_gradientn(limits = c(0,1),colours = c("#e5f5e0","#a1d99b","#31a354")) +
    ylab("") + xlab("") + coord_fixed() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)) +
    ggtitle(paste("Fraction of included SNVs ",p,sep=""))
  
  print(pp)
  
}

dev.off()
