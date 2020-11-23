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
library(ggrepel); packageVersion("ggrepel")

#theme
theme_set(theme_classic())

#load scripts
source('./scripts/plot_distances.R')
source('./scripts/color_palettes.R')

###################################
## defined functions
###################################
#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

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

BAC_pruned.ra.long$Order <- as.character(BAC_pruned.ra.long$Order)
BAC_pruned.ra.long$Order[is.na(BAC_pruned.ra.long$Order)] <- paste(BAC_pruned.ra.long$Class[is.na(BAC_pruned.ra.long$Order)],"uc", sep = "_")

BAC_pruned.ra.long$Family <- as.character(BAC_pruned.ra.long$Family)
BAC_pruned.ra.long$Family[is.na(BAC_pruned.ra.long$Family)] <- paste(BAC_pruned.ra.long$Order[is.na(BAC_pruned.ra.long$Family)],"uc", sep = "_")

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
BAC_pruned.ra.long %>% select(Group, layers, StationName, Type, Community, Class, Abundance)%>%
  group_by(Group, StationName, Type, Community,layers, Class) %>%
  summarize(class.abund= sum(Abundance)) -> tax_overview_by_stations


tax_overview_by_stations %>% group_by(Group, layers, Community, Class) %>%
        summarize(mean.abund= mean(class.abund),
                  se.abund = se(class.abund)) -> tax_overview_by_groups_layers

#overview of different famalies by layers
BAC_pruned.ra.long %>% select(Group, layers, StationName, Type, Community, Class, Order, Family, Abundance)%>%
  group_by(Group, StationName, Type, Community,layers, Class, Order,Family) %>%
  summarize(class.abund= sum(Abundance)) -> family_overview_by_stations


family_overview_by_stations %>% group_by(Group, layers, Community, Class, Order, Family) %>%
  summarize(mean.abund= mean(class.abund),
            se.abund = se(class.abund)) -> family_overview_by_groups_layers


#####################################
#PCA plot
#####################################
PS107_merged.prev.vst <- readRDS("./Data/PS107_merged_prev_vst.rds")

PS107.ord <- ordinate(PS107_merged.prev.vst, method = "RDA", distance = "eucledian")
PS107.ord.df <- plot_ordination(PS107_merged.prev.vst, PS107.ord, axes = c(1,2,3),justDF = TRUE)

#extract explained variance
PS107.ord.evals <- 100 * (PS107.ord$CA$eig/ sum(PS107.ord$CA$eig))
PS107.ord.df$ID <- rownames(PS107.ord.df)

PS107.ord.p <- ggplot(data = PS107.ord.df, aes(x =PC1, y=PC2, shape = Community))+
  geom_point(colour="black",size = 5)+
  geom_point(size = 4)+
  #stat_ellipse(data = PS107.ord.df, aes(x =PC1, y=PC2, group = Community), size = 1, type= "norm")+
  #geom_text(aes(label = StationName), colour = "black", nudge_y= -0.5,  size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(PS107.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(PS107.ord.evals[2], 2)), shape = "Depth", color = "Community")+
  scale_color_manual(values = c("Surface-10" = "lightblue", 
                                "Chl.max-20-30"="green",
                                "B.Chl.max-50"="darkgreen",
                                "Epipelagic-100"="blue",
                                "Mesopelagic-200"="darkblue",
                                "Mesopelagic-400"="darkblue")) +
  geom_label_repel(aes(label = paste(StationName,paste(Depth,"m", sep = "")),
                       fill = Group), color = 'white',
                   size = 3, seed = 1)+
  scale_fill_manual(values = c("in"= "red","out" = "blue"))+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")

ggsave("./figures/PCA_prev.pdf", PS107.ord.p, dpi = 300, 
       #width = 11.4, height = 23, 
       units = "cm")


#####################################
#PCA plot of upper 50 m
#####################################
PS107_upper50m <- subset_samples(PS107_merged.prev.vst, Type %in% c("Surface-10", 
                                                                    "Chl.max-20-30",
                                                                    "B.Chl.max-50"))

PS107_upper50m.ord <- ordinate(PS107_upper50m, method = "RDA", distance = "eucledian")
PS107_upper50m.ord.df <- plot_ordination(PS107_upper50m, PS107_upper50m.ord, axes = c(1,2,3),justDF = TRUE)

#extract explained variance
PS107_upper50m.ord.evals <- 100 * (PS107_upper50m.ord$CA$eig/ sum(PS107_upper50m.ord$CA$eig))
PS107_upper50m.ord.df$ID <- rownames(PS107_upper50m.ord.df)
PS107_upper50m.ord.df$Type <- factor(PS107_upper50m.ord.df$Type , levels=c("Surface-10", 
                                                                    "Chl.max-20-30",
                                                                    "B.Chl.max-50"))

PS107_upper50m.ord.p <- ggplot(data = PS107_upper50m.ord.df, aes(x =PC1, y=PC2, shape = Community, colour = Type))+
  geom_point(colour="black",size = 5)+
  geom_point(size = 4)+
  #geom_text_repel(aes(label = paste(StationName,paste(Depth,"m", sep = "")), colour = Group), fill = "white", size = 5, seed = 1, segment.colour = NA)+
  #stat_ellipse(data = PS107.ord.df, aes(x =PC1, y=PC2, group = Type), colour = "black",size = 1, type= "norm")+
  geom_text(aes(label = paste(StationName,paste(Depth,"m", sep = "")), colour = Group), nudge_y= -0.5,  size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(PS107.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(PS107.ord.evals[2], 2)), shape = "Depth", color = "Community")+
  scale_color_manual(values = c("Surface-10" = "lightblue", 
                                "Chl.max-20-30"="green",
                                "B.Chl.max-50"="darkgreen",
                                "in"= "red","out" = "blue")) +
  scale_fill_manual(values = c("in"= "red","out" = "blue"))+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")



ggsave("./figures/PCA_prev_surf.pdf", PS107_upper50m.ord.p, dpi = 300, 
       #width = 11.4, height = 23, 
       units = "cm")

#significance test
df <- as(sample_data(PS107_merged.prev.vst), "data.frame")
d <- phyloseq::distance(PS107_merged.prev.vst, "euclidean")
adonis_all <- adonis(d ~ Community + Group + Type, df)
adonis_all


#significance test on only surface samples
PS107_SRF <- subset_samples(PS107_merged.prev.vst, Type %in% c("Surface-10","Chl.max-20-30","B.Chl.max-50","Epipelagic-100"))

df <- as(sample_data(PS107_SRF), "data.frame")
d <- phyloseq::distance(PS107_SRF, "euclidean")
adonis_all <- adonis(d ~ Community + Group + Type, df)
adonis_all


#significance test below 50 m
PS107_deep <- subset_samples(PS107_merged.prev.vst, Type %in% c("Epipelagic-100","Mesopelagic-200","Mesopelagic-400"))

df <- as(sample_data(PS107_deep), "data.frame")
d <- phyloseq::distance(PS107_deep, "euclidean")
adonis_all <- adonis(d ~ Community + Group + Type, df)
adonis_all

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
