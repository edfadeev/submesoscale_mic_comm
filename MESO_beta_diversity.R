#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

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
PS107_merged.prev <- readRDS("./Data/PS107_merged_prev.rds")
PS107_merged.prev.vst <- readRDS("./Data/PS107_merged_prev_vst.rds")

#####################################
#Plot barplots of communities
#####################################
#transform data
BAC_pruned.ra <- transform_sample_counts(PS107_merged.prev, function(x) x / sum(x))
BAC_pruned.ra.long <- psmelt(BAC_pruned.ra)

#calculate abundance for each taxa
BAC_pruned.ra.long.agg <- aggregate(Abundance~StationName+Fraction+Type+Class, BAC_pruned.ra.long, FUN = "sum")
BAC_pruned.ra.long.agg$Abundance <- BAC_pruned.ra.long.agg$Abundance*100

#order of stations
levels(BAC_pruned.ra.long.agg$StationName) <- c("EG1","EG4","N5","N4","HG9", "HG4","HG2","HG1","S3")
#correct taxonomy
BAC_pruned.ra.long.agg$Class <- gsub("_unclassified","_uc",BAC_pruned.ra.long.agg$Class)

#remove below 2% ra
taxa_classes <- unique(BAC_pruned.ra.long.agg$Class)
BAC_pruned.ra.long.agg$Class[BAC_pruned.ra.long.agg$Abundance<2] <- "Other taxa"
BAC_pruned.ra.long.agg$Class <- factor(BAC_pruned.ra.long.agg$Class,
                                       levels=c(taxa_classes,"Other taxa"))

BAC_pruned.ra.long.agg.sub <- BAC_pruned.ra.long.agg[BAC_pruned.ra.long.agg$Type!="Sediment",]

#Plot 
barplots <- ggplot(BAC_pruned.ra.long.agg.sub, aes(x = StationName, y = Abundance, fill = Class)) + 
  facet_grid(Type~Fraction, space= "fixed") +
  geom_col()+
  #scale_fill_manual(values = phyla.col.PS107_shared) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme(legend.position="bottom")

ggsave("./figures/barplot_sed.pdf", barplots, dpi = 300, 
       width = 30, height = 30, 
       units = "cm")

#####################################
#PCoA plot
#####################################
#plot only water samples
PS107_merged.vst.no.sed <- subset_samples(PS107_merged.prev.vst,Type!="Sediment")

#calculate ordination
PS107.ord <- ordinate(PS107_merged.vst.no.sed, method = "RDA", distance = "eucledian")
PS107.ord.df <- plot_ordination(PS107_merged.vst.no.sed, PS107.ord, axes = c(1,2,3),justDF = TRUE)

#adjust grouping for clustering
PS107.ord.df$new_ordination <- paste(PS107.ord.df$Type, PS107.ord.df$Fraction, sep= ".")
PS107.ord.df$Region<- factor(
  PS107.ord.df$Region, 
  levels = c("EGC","WSC"))

#extract explained variance
PS107.ord.evals <- 100 * (PS107.ord$CA$eig/ sum(PS107.ord$CA$eig))
PS107.ord.df$ID <- rownames(PS107.ord.df)

PS107.ord.p <- ggplot(data = PS107.ord.df, aes(x =PC1, y =PC2, shape = Fraction, colour = Region, group = new_ordination))+
    geom_point(colour="black",size = 4)+
  geom_point(size = 3)+
  #geom_text(aes(label = Type), colour = "black", nudge_y= -1,  size=3)+
  labs(x = sprintf("PC1 [%s%%]", round(PS107.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(PS107.ord.evals[2], 2)), shape = "Fraction", color = "Origin")+
  stat_ellipse(colour = "black", size = 0.5)+
  scale_color_manual(values = c("EGC" = "blue", "WSC"="red")) +
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")

ggsave("./figures/Ordination_no_sed.pdf", PS107.ord.p, dpi = 300, 
    #width = 11.4, height = 23, 
     units = "cm")

#calculate dendograme
PS107_merged.dist <- dist(t(otu_table(PS107_merged.vst.no.sed)), method = "euclidean")

hc <- hclust(PS107_merged.dist, method = "ward.D")

plot(hc, labels = paste(sample_data(PS107_merged.vst.no.sed)$StationName,
                        sample_data(PS107_merged.vst.no.sed)$Fraction,
                        sample_data(PS107_merged.vst.no.sed)$Type))

#####################################
#Significance test
#####################################
df <- as(sample_data(PS107_merged.vst.no.sed), "data.frame")
d <- phyloseq::distance(PS107_merged.vst.no.sed, "euclidean")
adonis_all <- adonis(d ~ Fraction + Region + Type, df)
adonis_all


PS107_merged.vst.fl <- subset_samples(PS107_merged.vst.no.sed, Fraction == "FL") 
df.fl <- as(sample_data(PS107_merged.vst.fl), "data.frame")
d.fl <- phyloseq::distance(PS107_merged.vst.fl, "euclidean")
adonis_fl <- adonis(d.fl ~ Type + Region , df.fl)
adonis_fl

PS107_merged.vst.pa <- subset_samples(PS107_merged.vst.no.sed, Fraction == "PA") 
df.pa <- as(sample_data(PS107_merged.vst.pa), "data.frame")
d.pa <- phyloseq::distance(PS107_merged.vst.pa, "euclidean")
adonis_pa <- adonis(d.pa ~ Type + Region, df.pa)
adonis_pa

PS107_merged.vst.epi <- subset_samples(PS107_merged.vst.no.sed, Type == "SRF"| Type == "EPI")
PS107_merged.vst.epi <- subset_samples(PS107_merged.vst.epi, Fraction =="FL")
df <- as(sample_data(PS107_merged.vst.epi), "data.frame")
d <- phyloseq::distance(PS107_merged.vst.epi, "euclidean")
adonis_epi <- adonis(d ~ Region, df)
adonis_epi

PS107_merged.vst.epi <- subset_samples(PS107_merged.vst.no.sed, Type == "SRF"| Type == "EPI")
PS107_merged.vst.epi <- subset_samples(PS107_merged.vst.epi, Fraction =="PA")
df <- as(sample_data(PS107_merged.vst.epi), "data.frame")
d <- phyloseq::distance(PS107_merged.vst.epi, "euclidean")
adonis_epi <- adonis(d ~ Region, df)
adonis_epi

#check sediment association with regions
PS107_merged.vst.sed <- subset_samples(PS107_merged.prev.vst, Type == "Sediment")
df <- as(sample_data(PS107_merged.vst.sed), "data.frame")
d <- phyloseq::distance(PS107_merged.vst.sed, "euclidean")
adonis_sed <- adonis(d ~ Region, df)
adonis_sed

#####################################
#plot distances distribution between fractions
#####################################
sample_data(PS107_merged.vst.no.sed)$SampleID <- sample_names(PS107_merged.vst.no.sed)
wu.sd <- data.frame()

for (i in levels(droplevels(sample_data(PS107_merged.vst.no.sed)$Type))){
PS107_merged.vst.deep <- subset_samples(PS107_merged.vst.no.sed, Type == i)
PS107_merged.vst.deep <- prune_taxa(taxa_sums(PS107_merged.vst.deep)>0,PS107_merged.vst.deep)
# calc distances
wu = phyloseq::distance(PS107_merged.vst.deep, "euclidean")
wu.m = melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(PS107_merged.vst.deep) %>%
  select(SampleID, Fraction, Type, Source.Alt)

# combined distances with sample data
colnames(sd) = c("Var1", "Fraction", "Type", "Source.Alt")
test = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Fraction2", "Type2", "Source.Alt2")
test = left_join(test, sd, by = "Var2")

wu.sd <- rbind(wu.sd,test)
}

#filter out replication of distances
wu.sd_uni = wu.sd %>%
  filter(as.character(Fraction) != as.character(Fraction2)) %>%
  filter(Fraction =="PA")

wu.sd_uni$Type <- factor(wu.sd_uni$Type, levels = rev(levels(wu.sd_uni$Type)))

#plot
FL_PA_dist.p <- ggplot(wu.sd_uni, aes(x = Type, y = value)) +
  theme_classic(base_size = 12) +
  labs(x = "Water layer")+
  geom_boxplot(outlier.color = NULL, notch = FALSE)+
  #geom_jitter(position=position_jitter(0.2), alpha =0.2, colour = "black")+
  scale_fill_manual(values = c("SRF" = "yellow", "EPI"= "pink", "MESO"="darkblue","BATHY"= "gray")) +
  #ggtitle(paste0("Distance Metric = ", "euclidean"))+
  coord_flip()+
  geom_signif(comparisons = list(c("SRF", "EPI"),c("MESO","BATHY"),c("EPI", "MESO")), 
              map_signif_level=TRUE, test = "wilcox.test")+
  theme(legend.position = "none")

#####################################
#plot distances distribution between depths
#####################################
sample_data(PS107_merged.vst.no.sed)$SampleID <- sample_names(PS107_merged.vst.no.sed)
wu.sd_frac <- data.frame()

for (i in levels(droplevels(sample_data(PS107_merged.vst.no.sed)$Fraction))){
  
  PS107_merged.vst.deep <- subset_samples(PS107_merged.vst.no.sed, Fraction == i)
  PS107_merged.vst.deep <- prune_taxa(taxa_sums(PS107_merged.vst.deep)>0,PS107_merged.vst.deep)
  # calc distances
  wu = phyloseq::distance(PS107_merged.vst.deep, "euclidean")
  wu.m = melt(as.matrix(wu))
  
  # remove self-comparisons
  wu.m = wu.m %>%
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor,as.character)
  
  # get sample data (S4 error OK and expected)
  sd = sample_data(PS107_merged.vst.deep) %>%
    select(SampleID, Fraction, Fraction, Type)
  
  # combined distances with sample data
  colnames(sd) = c("Var1", "Fraction", "Type")
  test = left_join(wu.m, sd, by = "Var1")
  
  colnames(sd) = c("Var2", "Fraction2", "Type2")
  test = left_join(test, sd, by = "Var2")
  
  wu.sd_frac <- rbind(wu.sd_frac,test)
}

#filter out replication of distances
wu.sd_frac.uni = wu.sd_frac %>%
  filter(as.character(Type) != as.character(Type2)) 


frac_dist.p <- ggplot(wu.sd_frac.uni, aes(x = Fraction, y = value)) +
  theme_classic(base_size = 12) +
  labs(y = "Euclidean distance")+
  geom_boxplot(outlier.color = NULL, notch = FALSE)+
  #geom_jitter(position=position_jitter(0.2), alpha= 0.2, colour = "black")+
  #scale_fill_manual(values = c("DCM" = "yellow", "EPI"= "pink", "MESO"="darkblue","BATHY"= "gray")) +
  #ggtitle(paste0("Distance Metric = ", "euclidean"))+
  geom_signif(comparisons = list(c("FL", "PA")), 
              map_signif_level=TRUE, test = "wilcox.test")+
  theme(legend.position = "none")



#combined plot
plot_grid(FL_PA_dist.p, frac_dist.p, labels = c("A","B"), ncol = 2, align = "h")

ggsave("./figures/Figure-beta-no-leg.pdf", 
       plot = ggplot2::last_plot(),
       scale = 1,
       units = "cm",
       width = 17.8,
       #height = 17.4,
       dpi = 300)



#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
