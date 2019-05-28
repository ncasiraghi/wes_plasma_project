wd <- "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_mutation_load_tables"
setwd(wd)

MutationLoadTable <- "/elaborazioni/sharedCO/Abemus_data_analysis/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/Mutect_mutation_load_tables/MutationLoadTable_Mutect_afn_0.01_aft_0.05_mincov_30_mintc_0.1.RData"

# load the MutationLoadTable
load(MutationLoadTable)
final$ratio_median_all <- NA

final$ratio_median_all[which(final$Type=="Tissue")] <- final$n.snvs[which(final$Type=="Tissue")]/median(final$n.snvs[which(final$Type=="Tissue")])

final$ratio_median_all[which(final$Type=="Plasma")] <- final$n.snvs[which(final$Type=="Plasma")]/median(final$n.snvs[which(final$Type=="Plasma")])

save(final,file = "ratio_mutload_medianload_Mutect_afn_0.01_aft_0.05_mincov_30_mintc_0.1.RData",compress = T)