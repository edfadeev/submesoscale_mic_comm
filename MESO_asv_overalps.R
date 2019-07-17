#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
#load libraries
library("VennDiagram"); packageVersion("VennDiagram")
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("olsrr"); packageVersion("olsrr")
library("cowplot"); packageVersion("cowplot")
library("venn"); packageVersion("venn")

#load colour palettes
source('./Scripts/color_palettes.R')

#####################################
#Load phyloseq object
####################################
PS107_merged <-  readRDS("./Data/PS107_merged.rds")

#####################################
#Calculate OTU overlapos for each depth
####################################

for (frac in sample_data(PS107_merged)$Community){
  PS107_merged.frac<- subset_samples(PS107_merged, Community == frac)
    for (type in sample_data(PS107_merged)$Type){
      PS107_merged.type <- subset_samples(PS107_merged.frac, Type== type)
      PS107_merged.type <- prune_taxa(taxa_sums(PS107_merged.type)>0,PS107_merged.type)
      y <- list()
        for (station in sample_data(PS107_merged.type)$StationName){
          T1_comm<- subset_samples(PS107_merged.type, StationName == station)
          T1_comm <- prune_taxa(taxa_sums(T1_comm)>0,T1_comm)
          y[[station]] <- taxa_names(T1_comm)
        }
      assign(paste(type, frac, sep = "_"), y)
    }
}


pdf("./figures/ASVs_overlap_FL.pdf")
venn(`Surface-10_FL`, snames = names(`Surface-10_FL`), ilab=TRUE, zcolor = "style")
venn(`Chl.max-20-30_FL`, snames = names(`Chl.max-20-30_FL`), ilab=TRUE, zcolor = "style")
venn(`B.Chl.max-50_FL`, snames = names(`B.Chl.max-50_FL`), ilab=TRUE, zcolor = "style")
venn(`Epipelagic-100_FL`, snames = names(`Epipelagic-100_FL`), ilab=TRUE, zcolor = "style")
venn(`Mesopelagic-200_FL`, snames = names(`Mesopelagic-200_FL`), ilab=TRUE, zcolor = "style")
venn(`Mesopelagic-400_FL`, snames = names(`Mesopelagic-400_FL`), ilab=TRUE, zcolor = "style")
dev.off()

pdf("./figures/ASVs_overlap_PA.pdf")
venn(`Surface-10_PA`, snames = names(`Surface-10_PA`), ilab=TRUE, zcolor = "style")
venn(`Chl.max-20-30_PA`, snames = names(`Chl.max-20-30_PA`), ilab=TRUE, zcolor = "style")
venn(`B.Chl.max-50_PA`, snames = names(`B.Chl.max-50_PA`), ilab=TRUE, zcolor = "style")
venn(`Epipelagic-100_PA`, snames = names(`Epipelagic-100_PA`), ilab=TRUE, zcolor = "style")
venn(`Mesopelagic-200_PA`, snames = names(`Mesopelagic-200_PA`), ilab=TRUE, zcolor = "style")
venn(`Mesopelagic-400_PA`, snames = names(`Mesopelagic-400_PA`), ilab=TRUE, zcolor = "style")
dev.off()
