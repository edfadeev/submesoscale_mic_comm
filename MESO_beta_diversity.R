#set working directory in Windows
# wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(wd)

#load libraries
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
library(reshape2); packageVersion("reshape2")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(phyloseq); packageVersion("phyloseq")
library(ggsignif); packageVersion("ggsignif")

#theme
theme_set(theme_classic())

#load scripts
source('./scripts/plot_distances.R')
source('./scripts/color_palettes.R')

#Load datatset
PS107_merged <- readRDS("./Data/PS107_merged.rds")
PS107_merged.prev <- readRDS("./Data/PS107_merged_prev.rds")


#####################################
#Plot barplots of communities
#####################################
#transform data
BAC_pruned.ra <- transform_sample_counts(PS107_merged.prev, function(x) x / sum(x))
BAC_pruned.ra.long <- psmelt(BAC_pruned.ra)
BAC_pruned.ra.long$Abundance <- BAC_pruned.ra.long$Abundance*100

#fix unclassified lineages 
BAC_pruned.ra.long$Class <- as.character(BAC_pruned.ra.long$Class)
BAC_pruned.ra.long$Class[is.na(BAC_pruned.ra.long$Class)] <- paste(BAC_pruned.ra.long$Phylum[is.na(BAC_pruned.ra.long$Class)],"uc", sep = "_")

#calculate abundance for each Class
BAC_pruned.ra.long %>% select(StationName,Community,Type,Class,Abundance)%>%
                        group_by(StationName,Community,Type,Class) %>%
                        summarize(Abund.total= sum(Abundance)) -> BAC_pruned.ra.long.agg

#order of stations
BAC_pruned.ra.long.agg$StationName <- factor(BAC_pruned.ra.long.agg$StationName, 
                                             levels = c("T4","T1","T2","T5","T3"))
BAC_pruned.ra.long.agg$Type <- factor(BAC_pruned.ra.long.agg$Type,
                                      levels = c("Surface-10","Chl.max-20-30","B.Chl.max-50","Epipelagic-100","Mesopelagic-200","Mesopelagic-400"))

#remove below 1% ra
taxa_classes <- unique(BAC_pruned.ra.long.agg$Class)
BAC_pruned.ra.long.agg$Class[BAC_pruned.ra.long.agg$Abund.total<1] <- "Other taxa"
BAC_pruned.ra.long.agg$Class <- factor(BAC_pruned.ra.long.agg$Class,
                                       levels=c(taxa_classes,"Other taxa"))

#Plot 
barplots <- ggplot(BAC_pruned.ra.long.agg, aes(x = StationName, y = Abund.total, fill = Class)) + 
  facet_grid(Type~Community, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = sample(tol21rainbow)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme(legend.position="bottom")

ggsave("./figures/barplot-prev.pdf", barplots, dpi = 300, 
       width = 30, height = 30, 
       units = "cm")

#overview of different groups by layers
BAC_pruned.ra.long %>% select(Group, layers, StationName, Type, Community,Class, Abundance)%>%
  group_by(Group, StationName, Type, Community,layers, Class) %>%
  summarize(class.abund= sum(Abundance)) -> tax_overview_by_stations


tax_overview_by_stations %>% group_by(Group, layers, Community, Class) %>%
        summarize(mean.abund= mean(class.abund),
                  se.abund = se(class.abund)) -> tax_overview_by_groups_layers


#####################################
#PCA plot
#####################################
PS107.ord <- ordinate(PS107_merged.prev.vst, method = "RDA")
PS107.ord.df <- plot_ordination(PS107_merged.prev.vst, PS107.ord, axes = c(1,2,3),justDF = TRUE)

#adjust grouping for clustering
PS107.ord.df$new_ordination <- paste(PS107.ord.df$StationName, PS107.ord.df$Depth,PS107.ord.df$Community, sep= ".")

#extract explained variance
PS107.ord.evals <- 100 * (PS107.ord$CA$eig/ sum(PS107.ord$CA$eig))
PS107.ord.df$ID <- rownames(PS107.ord.df)

PS107.ord.p <- ggplot(data = PS107.ord.df, aes(x =PC1, y=PC2, shape = Type, colour = Community))+
  geom_point(colour="black",size = 4)+
  geom_point(size = 3)+
  #geom_polygon(data=PS107.ord.df,aes(x=NMDS1,y=NMDS2,fill=Type,group=Type),alpha=0.30) +
  geom_text(aes(label = StationName), colour = "black", nudge_y= -0.5,  size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(PS107.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(PS107.ord.evals[2], 2)), shape = "Depth", color = "Community")+
  #annotate(geom="text", size = 5, x=-1.2, y=0.7, label= paste("Stress =", round(PS107.ord$stress, 3), sep = " "),color="black")+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")

ggsave("./figures/PCA_prev.pdf", PS107.ord.p, dpi = 300, 
       #width = 11.4, height = 23, 
       units = "cm")


#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
