# #set working directory
# wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(wd)

#load libraries
library("phyloseq"); packageVersion("phyloseq")
library("venn"); packageVersion("venn")
library("UpSetR"); packageVersion("UpSetR")


library("VennDiagram"); packageVersion("VennDiagram")

library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("olsrr"); packageVersion("olsrr")
library("cowplot"); packageVersion("cowplot")



#load colour palettes
source('./Scripts/color_palettes.R')
source("./Scripts/pres_abs_matrix.R")
#####################################
#Load phyloseq object
####################################
PS107_merged.prev <-  readRDS("./Data/PS107_merged_prev.rds")

#####################################
#venn Diagram of bacterial ASV overlap between the different groups in upper 50 m
#####################################
#FL
BAC_FL <- subset_samples(PS107_merged.prev, Community == "FL")
BAC_FL <- prune_taxa(taxa_sums(BAC_FL)>0,BAC_FL)

#up
BAC_FL.up <- subset_samples(BAC_FL, layers == "up")
BAC_FL.up <- prune_taxa(taxa_sums(BAC_FL.up)>0,BAC_FL.up)

BAC_FL.up.T1 <- subset_samples(BAC_FL.up, StationName == "T1")
BAC_FL.up.T1 <- prune_taxa(taxa_sums(BAC_FL.up.T1)>0,BAC_FL.up.T1)

BAC_FL.up.T2 <- subset_samples(BAC_FL.up, StationName == "T2")
BAC_FL.up.T2 <- prune_taxa(taxa_sums(BAC_FL.up.T2)>0,BAC_FL.up.T2)

BAC_FL.up.T3 <- subset_samples(BAC_FL.up, StationName == "T3")
BAC_FL.up.T3 <- prune_taxa(taxa_sums(BAC_FL.up.T3)>0,BAC_FL.up.T3)

BAC_FL.up.T4 <- subset_samples(BAC_FL.up, StationName == "T4")
BAC_FL.up.T4 <- prune_taxa(taxa_sums(BAC_FL.up.T4)>0,BAC_FL.up.T4)

BAC_FL.up.T5 <- subset_samples(BAC_FL.up, StationName == "T5")
BAC_FL.up.T5 <- prune_taxa(taxa_sums(BAC_FL.up.T5)>0,BAC_FL.up.T5)

#make a list
y <- list()
y[["FL-T1-UP-IN"]] <- as.character(row.names(otu_table(BAC_FL.up.T1)))
y[["FL-T2-UP-IN"]] <- as.character(row.names(otu_table(BAC_FL.up.T2)))
y[["FL-T3-UP-OUT"]] <- as.character(row.names(otu_table(BAC_FL.up.T3)))
y[["FL-T4-UP-OUT"]] <- as.character(row.names(otu_table(BAC_FL.up.T4)))
y[["FL-T5-UP-IN"]] <- as.character(row.names(otu_table(BAC_FL.up.T5)))


#PA
BAC_PA <- subset_samples(PS107_merged.prev, Community == "PA")
BAC_PA <- prune_taxa(taxa_sums(BAC_PA)>0,BAC_PA)

#up
BAC_PA.up <- subset_samples(BAC_PA, layers == "up")
BAC_PA.up <- prune_taxa(taxa_sums(BAC_PA.up)>0,BAC_PA.up)

BAC_PA.up.T1 <- subset_samples(BAC_PA.up, StationName == "T1")
BAC_PA.up.T1 <- prune_taxa(taxa_sums(BAC_PA.up.T1)>0,BAC_PA.up.T1)

BAC_PA.up.T2 <- subset_samples(BAC_PA.up, StationName == "T2")
BAC_PA.up.T2 <- prune_taxa(taxa_sums(BAC_PA.up.T2)>0,BAC_PA.up.T2)

BAC_PA.up.T3 <- subset_samples(BAC_PA.up, StationName == "T3")
BAC_PA.up.T3 <- prune_taxa(taxa_sums(BAC_PA.up.T3)>0,BAC_PA.up.T3)

BAC_PA.up.T4 <- subset_samples(BAC_PA.up, StationName == "T4")
BAC_PA.up.T4 <- prune_taxa(taxa_sums(BAC_PA.up.T4)>0,BAC_PA.up.T4)

BAC_PA.up.T5 <- subset_samples(BAC_PA.up, StationName == "T5")
BAC_PA.up.T5 <- prune_taxa(taxa_sums(BAC_PA.up.T5)>0,BAC_PA.up.T5)

#make a list
z <- list()
z[["PA-T1-UP-IN"]] <- as.character(row.names(otu_table(BAC_PA.up.T1)))
z[["PA-T2-UP-IN"]] <- as.character(row.names(otu_table(BAC_PA.up.T2)))
z[["PA-T3-UP-OUT"]] <- as.character(row.names(otu_table(BAC_PA.up.T3)))
z[["PA-T4-UP-OUT"]] <- as.character(row.names(otu_table(BAC_PA.up.T4)))
z[["PA-T5-UP-IN"]] <- as.character(row.names(otu_table(BAC_PA.up.T5)))


#plot
venn(y, snames = names(y), ilab=TRUE, zcolor = "style")

venn(z, snames = names(y), ilab=TRUE, zcolor = "style")


#####################################
#upset diagrame for the entire water column
#####################################
#FL
#down
BAC_FL.down <- subset_samples(BAC_FL, layers == "down")
BAC_FL.down <- prune_taxa(taxa_sums(BAC_FL.down)>0,BAC_FL.down)

BAC_FL.down.T1 <- subset_samples(BAC_FL.down, StationName == "T1")
BAC_FL.down.T1 <- prune_taxa(taxa_sums(BAC_FL.down.T1)>0,BAC_FL.down.T1)

BAC_FL.down.T2 <- subset_samples(BAC_FL.down, StationName == "T2")
BAC_FL.down.T2 <- prune_taxa(taxa_sums(BAC_FL.down.T2)>0,BAC_FL.down.T2)

BAC_FL.down.T3 <- subset_samples(BAC_FL.down, StationName == "T3")
BAC_FL.down.T3 <- prune_taxa(taxa_sums(BAC_FL.down.T3)>0,BAC_FL.down.T3)

BAC_FL.down.T4 <- subset_samples(BAC_FL.down, StationName == "T4")
BAC_FL.down.T4 <- prune_taxa(taxa_sums(BAC_FL.down.T4)>0,BAC_FL.down.T4)

BAC_FL.down.T5 <- subset_samples(BAC_FL.down, StationName == "T5")
BAC_FL.down.T5 <- prune_taxa(taxa_sums(BAC_FL.down.T5)>0,BAC_FL.down.T5)

y[["FL-T1-down-IN"]] <- as.character(row.names(otu_table(BAC_FL.down.T1)))
y[["FL-T2-down-IN"]] <- as.character(row.names(otu_table(BAC_FL.down.T2)))
y[["FL-T3-down-OUT"]] <- as.character(row.names(otu_table(BAC_FL.down.T3)))
y[["FL-T4-down-OUT"]] <- as.character(row.names(otu_table(BAC_FL.down.T4)))
y[["FL-T5-down-IN"]] <- as.character(row.names(otu_table(BAC_FL.down.T5)))

#generate overlap matrix
otu_overlaps <- pres_abs_matrix(y)    
otu_overlaps$OTU <- rownames(otu_overlaps)

#####################################
#Explore overlaping OTU
#####################################
taxonomy <- as.data.frame(tax_table(PS107_merged.prev))
taxonomy$OTU <- rownames(taxonomy)

otu_overlaps_merged <- full_join(taxonomy,otu_overlaps, by = c("OTU"))

#set metadata
sets <- names(y)

metadata <- as.data.frame(cbind(sets, c(rep("FL",10)),
                                rep(c("T1","T2","T3","T4","T5"),2),
                                rep(c(rep("up",5),rep("down",5)))))
names(metadata) <- c("sets", "Fraction","Station","Depth")


upset(otu_overlaps_merged, number.angles = 30,
      sets = as.vector(metadata$sets),
      #keep.order = TRUE, 
      #mainbar.y.max = 900,
      mainbar.y.label = "No. of overlaping OTU",
      order.by = "freq")#empty.intersections = "on")

#PA
#down
BAC_PA.down <- subset_samples(BAC_PA, layers == "down")
BAC_PA.down <- prune_taxa(taxa_sums(BAC_PA.down)>0,BAC_PA.down)

BAC_PA.down.T1 <- subset_samples(BAC_PA.down, StationName == "T1")
BAC_PA.down.T1 <- prune_taxa(taxa_sums(BAC_PA.down.T1)>0,BAC_PA.down.T1)

BAC_PA.down.T2 <- subset_samples(BAC_PA.down, StationName == "T2")
BAC_PA.down.T2 <- prune_taxa(taxa_sums(BAC_PA.down.T2)>0,BAC_PA.down.T2)

BAC_PA.down.T3 <- subset_samples(BAC_PA.down, StationName == "T3")
BAC_PA.down.T3 <- prune_taxa(taxa_sums(BAC_PA.down.T3)>0,BAC_PA.down.T3)

BAC_PA.down.T4 <- subset_samples(BAC_PA.down, StationName == "T4")
BAC_PA.down.T4 <- prune_taxa(taxa_sums(BAC_PA.down.T4)>0,BAC_PA.down.T4)

BAC_PA.down.T5 <- subset_samples(BAC_PA.down, StationName == "T5")
BAC_PA.down.T5 <- prune_taxa(taxa_sums(BAC_PA.down.T5)>0,BAC_PA.down.T5)

z[["PA-T1-down-IN"]] <- as.character(row.names(otu_table(BAC_PA.down.T1)))
z[["PA-T2-down-IN"]] <- as.character(row.names(otu_table(BAC_PA.down.T2)))
z[["PA-T3-down-OUT"]] <- as.character(row.names(otu_table(BAC_PA.down.T3)))
z[["PA-T4-down-OUT"]] <- as.character(row.names(otu_table(BAC_PA.down.T4)))
z[["PA-T5-down-IN"]] <- as.character(row.names(otu_table(BAC_PA.down.T5)))


#generate overlap matrix
otu_overlaps <- pres_abs_matrix(z)    
otu_overlaps$OTU <- rownames(otu_overlaps)

#####################################
#Explore overlaping OTU
#####################################
taxonomy <- as.data.frame(tax_table(PS107_merged.prev))
taxonomy$OTU <- rownames(taxonomy)

otu_overlaps_merged <- full_join(taxonomy,otu_overlaps, by = c("OTU"))

#set metadata
sets <- names(z)

metadata <- as.data.frame(cbind(sets, c(rep("PA",10)),
                                rep(c("T1","T2","T3","T4","T5"),2),
                                rep(c(rep("up",5),rep("down",5)))))
names(metadata) <- c("sets", "Fraction","Station","Depth")


upset(otu_overlaps_merged, number.angles = 30,
      sets = as.vector(metadata$sets),
      #keep.order = TRUE, 
      #mainbar.y.max = 900,
      mainbar.y.label = "No. of overlaping OTU",
      order.by = "freq")#empty.intersections = "on")

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
